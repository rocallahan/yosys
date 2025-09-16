/*
 *  yosys -- Yosys Open SYnthesis Suite
 *
 *  Copyright (C) 2012  Claire Xenia Wolf <claire@yosyshq.com>
 *
 *  Permission to use, copy, modify, and/or distribute this software for any
 *  purpose with or without fee is hereby granted, provided that the above
 *  copyright notice and this permission notice appear in all copies.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 *  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 *  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 *  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 *  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 *  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 *  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

#include "kernel/register.h"
#include "kernel/ffinit.h"
#include "kernel/sigtools.h"
#include "kernel/log.h"
#include "kernel/celltypes.h"
#include "kernel/threading.h"
#include "kernel/trtlil.h"
#include "libs/sha1/sha1.h"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <array>


USING_YOSYS_NAMESPACE
PRIVATE_NAMESPACE_BEGIN

template <typename T, typename U>
inline Hasher hash_pair(const T &t, const U &u) { return hash_ops<std::pair<T, U>>::hash(t, u); }

struct ComputeCellHashes
{
	int cell_index_begin;
	int cell_index_end;
};

struct CellHash
{
	int cell_index;
	Hasher::hash_t hash_value;
};

struct ComputeCellHashesOut
{
	// Entry i contains the hashes where hash_value % bucketed_cell_hashes.size() == i
	std::vector<std::vector<CellHash>> bucketed_cell_hashes;
};

struct FindDuplicateCells
{
	std::vector<std::vector<std::vector<CellHash>>> &bucketed_cell_hashes;
};

struct DuplicateCell
{
	int remove_cell;
	int keep_cell;
};

struct FindDuplicateCellsOut
{
	std::vector<DuplicateCell> duplicates;
};

struct OptMergeThreadWorker
{
	TRTLIL::RModule module;
	const TRTLIL::TSigMap &assign_map;
	TRTLIL::RFfInitVals init_vals;
	TRTLIL::RCellTypes cell_types;
	int workers;
	bool mode_share_all;
	bool mode_keepdc;

	static Hasher hash_pmux_in(const TRTLIL::TSigSpec& sig_s, const TRTLIL::TSigSpec& sig_b, Hasher h)
	{
		int s_width = GetSize(sig_s);
		int width = GetSize(sig_b) / s_width;

		hashlib::commutative_hash comm;

		auto sig_s_it = sig_s.begin();
		for (int i = 0; i < s_width; ++i, ++sig_s_it)
			comm.eat(hash_pair(*sig_s_it, sig_b.extract(i*width, width)));

		return comm.hash_into(h);
	}

	static void sort_pmux_conn(dict<TRTLIL::RIdString, TRTLIL::TSigSpec> &conn)
	{
		const TRTLIL::TSigSpec &sig_s = conn.at(ID::S);
		const TRTLIL::TSigSpec &sig_b = conn.at(ID::B);

		int s_width = GetSize(sig_s);
		int width = GetSize(sig_b) / s_width;

		vector<pair<TRTLIL::TSigBit, TRTLIL::TSigSpec>> sb_pairs;
		auto s_it = sig_s.begin();
		for (int i = 0; i < s_width; ++i, ++s_it)
			sb_pairs.push_back({*s_it, sig_b.extract(i*width, width)});

		std::sort(sb_pairs.begin(), sb_pairs.end());

		conn[ID::S] = TRTLIL::TSigSpec();
		conn[ID::B] = TRTLIL::TSigSpec();

		for (auto &it : sb_pairs) {
			conn[ID::S].append(it.first);
			conn[ID::B].append(it.second);
		}
	}

	Hasher hash_cell_inputs(TRTLIL::RCell cell, Hasher h) const
	{
		// TODO: when implemented, use celltypes to match:
		// (builtin || stdcell) && (unary || binary) && symmetrical
		if (cell.get_type().in(ID($and), ID($or), ID($xor), ID($xnor), ID($add), ID($mul),
				ID($logic_and), ID($logic_or), ID($_AND_), ID($_OR_), ID($_XOR_))) {
			hashlib::commutative_hash comm;
			comm.eat(cell.get_port_mapped(ID::A, assign_map));
			comm.eat(cell.get_port_mapped(ID::B, assign_map));
			h = comm.hash_into(h);
		} else if (cell.get_type().in(ID($reduce_xor), ID($reduce_xnor))) {
			TRTLIL::TSigSpec a = cell.get_port_mapped(ID::A, assign_map);
			a.sort();
			h = a.hash_into(h);
		} else if (cell.get_type().in(ID($reduce_and), ID($reduce_or), ID($reduce_bool))) {
			TRTLIL::TSigSpec a = cell.get_port_mapped(ID::A, assign_map);
			a.sort_and_unify();
			h = a.hash_into(h);
		} else if (cell.get_type() == ID($pmux)) {
			TRTLIL::TSigSpec sig_s = cell.get_port_mapped(ID::S, assign_map);
			TRTLIL::TSigSpec sig_b = cell.get_port_mapped(ID::B, assign_map);
			h = hash_pmux_in(sig_s, sig_b, h);
			h = cell.get_port_mapped(ID::A, assign_map).hash_into(h);
		} else {
			hashlib::commutative_hash comm;
			for (const auto& [port, sig] : cell.mapped_ports(assign_map)) {
				if (cell.output(port))
					continue;
				comm.eat(hash_pair(port, sig));
			}
			h = comm.hash_into(h);
			if (cell.is_builtin_ff())
				h = cell.get_port_init_vals(ID::Q, init_vals).hash_into(h);
		}
		return h;

	}

	static Hasher hash_cell_parameters(TRTLIL::RCell cell, Hasher h)
	{
		hashlib::commutative_hash comm;
		for (const auto& param : cell.parameters()) {
			comm.eat(param);
		}
		return comm.hash_into(h);
	}

	Hasher hash_cell_function(TRTLIL::RCell cell, Hasher h) const
	{
		h.eat(cell.get_type());
		h = hash_cell_inputs(cell, h);
		h = hash_cell_parameters(cell, h);
		return h;
	}

	bool compare_cell_parameters_and_connections(TRTLIL::RCell cell1, TRTLIL::RCell cell2) const
	{
		if (cell1.same_cell(cell2)) return true;
		if (cell1.get_type() != cell2.get_type()) return false;

		if (cell1.parameters() != cell2.parameters())
			return false;

		TRTLIL::RCell::MappedPorts cell1_ports = cell1.mapped_ports(assign_map);
		int cell2_ports_size = cell2.mapped_ports(assign_map).size();
		if (cell1_ports.size() != cell2_ports_size)
			return false;
		for (const auto &it : cell1_ports)
			if (!cell2.has_port(it.first))
				return false;

		dict<TRTLIL::RIdString, TRTLIL::TSigSpec> conn1, conn2;
		conn1.reserve(cell1_ports.size());
		conn2.reserve(cell2_ports_size);

		for (const auto &it : cell1_ports) {
			if (cell1.output(it.first)) {
				if (it.first == ID::Q && cell1.is_builtin_ff()) {
					// For the 'Q' output of state elements,
					//   use the (* init *) attribute value
					conn1[it.first] = TRTLIL::TSigSpec::build(cell1.get_port_init_vals(ID::Q, init_vals));
					conn2[it.first] = TRTLIL::TSigSpec::build(cell2.get_port_init_vals(ID::Q, init_vals));
				}
				else {
					conn1[it.first] = TRTLIL::TSigSpec();
					conn2[it.first] = TRTLIL::TSigSpec();
				}
			}
			else {
				conn1[it.first] = it.second;
				conn2[it.first] = cell2.get_port_mapped(it.first, assign_map);
			}
		}

		if (cell1.get_type().in(ID($and), ID($or), ID($xor), ID($xnor), ID($add), ID($mul),
				ID($logic_and), ID($logic_or), ID($_AND_), ID($_OR_), ID($_XOR_))) {
			if (conn1.at(ID::A) < conn1.at(ID::B)) {
				std::swap(conn1[ID::A], conn1[ID::B]);
			}
			if (conn2.at(ID::A) < conn2.at(ID::B)) {


				std::swap(conn2[ID::A], conn2[ID::B]);
			}
		} else
		if (cell1.get_type().in(ID($reduce_xor), ID($reduce_xnor))) {
			conn1[ID::A].sort();
			conn2[ID::A].sort();
		} else
		if (cell1.get_type().in(ID($reduce_and), ID($reduce_or), ID($reduce_bool))) {
			conn1[ID::A].sort_and_unify();
			conn2[ID::A].sort_and_unify();
		} else
		if (cell1.get_type() == ID($pmux)) {
			sort_pmux_conn(conn1);
			sort_pmux_conn(conn2);
		}

		return conn1 == conn2;
	}

	bool has_dont_care_initval(TRTLIL::RCell cell) const
	{
		if (!cell.is_builtin_ff())
			return false;

		return !cell.get_port_init_vals(ID::Q, init_vals).is_fully_def();
	}

	OptMergeThreadWorker(const RTLIL::Module &module, const GenericFfInitVals<TRTLIL::TSigMap> &init_vals,
			const TRTLIL::TSigMap &assign_map, const CellTypes &cell_types, int workers,
			bool mode_share_all, bool mode_keepdc) :
		module(module), assign_map(assign_map), init_vals(init_vals), cell_types(cell_types),
		workers(workers), mode_share_all(mode_share_all), mode_keepdc(mode_keepdc)
	{
	}

	ComputeCellHashesOut compute_cell_hashes(const ComputeCellHashes &in) const
	{
		std::vector<std::vector<CellHash>> bucketed_cell_hashes(workers);
		for (int cell_index = in.cell_index_begin; cell_index < in.cell_index_end; ++cell_index) {
			TRTLIL::RCell cell = module.cell_at(cell_index);
			if (!module.selected(cell))
				continue;
			if (cell.get_type().in(ID($meminit), ID($meminit_v2), ID($mem), ID($mem_v2))) {
				// Ignore those for performance: meminit can have an excessively large port,
				// mem can have an excessively large parameter holding the init data
				continue;
			}
			if (cell.get_type() == ID($scopeinfo))
				continue;
			if (mode_keepdc && has_dont_care_initval(cell))
				continue;
			if (!cell.known())
				continue;
			if (!mode_share_all && !cell_types.cell_known(cell.get_type()))
				continue;

			Hasher::hash_t h = hash_cell_function(cell, Hasher()).yield();
			int bucket_index = h % workers;
			bucketed_cell_hashes[bucket_index].push_back({cell_index, h});
		}
		return {std::move(bucketed_cell_hashes)};
	}

	FindDuplicateCellsOut find_duplicate_cells(int index, const FindDuplicateCells &in) const
	{
		// We keep a set of known cells. They're hashed with our hash_cell_function
		// and compared with our compare_cell_parameters_and_connections.
		// Both need to capture OptMergeThreadWorker.
		struct CellHashOp {
			std::size_t operator()(const CellHash &c) const {
				return (std::size_t)c.hash_value;
			}
		};
		struct CellEqualOp {
			const OptMergeThreadWorker& worker;
			CellEqualOp(const OptMergeThreadWorker& w) : worker(w) {}
			bool operator()(const CellHash &lhs, const CellHash &rhs) const {
				return worker.compare_cell_parameters_and_connections(
						worker.module.cell_at(lhs.cell_index),
						worker.module.cell_at(rhs.cell_index));
			}
		};
		std::unordered_set<
			CellHash,
			CellHashOp,
			CellEqualOp> known_cells(0, CellHashOp(), CellEqualOp(*this));

		std::vector<DuplicateCell> duplicates;
		for (const std::vector<std::vector<CellHash>> &buckets : in.bucketed_cell_hashes) {
			// Clear out our buckets as we go. This keeps the work of deallocation
			// off the main thread.
			std::vector<CellHash> bucket = std::move(buckets[index]);
			for (CellHash c : bucket) {
				auto [cell_in_map, inserted] = known_cells.insert(c);
				if (inserted)
					continue;
				CellHash map_c = *cell_in_map;
				if (module.cell_at(c.cell_index).has_keep_attr()) {
					if (module.cell_at(map_c.cell_index).has_keep_attr())
						continue;
					known_cells.erase(map_c);
					known_cells.insert(c);
					std::swap(c, map_c);
				}
				duplicates.push_back({c.cell_index, map_c.cell_index});
			}
		}
		return {duplicates};
	}
};

template <typename T>
void initialize_queues(std::vector<ConcurrentQueue<T>> &queues, int size) {
	queues.reserve(size);
	for (int i = 0; i < size; ++i)
		queues.emplace_back(1);
}

struct OptMergeWorker
{
	int total_count;

	OptMergeWorker(RTLIL::Module *module, bool mode_nomux, bool mode_share_all, bool mode_keepdc) :
		total_count(0)
	{
		TRTLIL::TSigMap assign_map = TRTLIL::TSigMap::build(*module);
		GenericFfInitVals<TRTLIL::TSigMap> init_vals;
		init_vals.set(&assign_map, module);

		CellTypes ct;
		ct.setup_internals();
		ct.setup_internals_mem();
		ct.setup_stdcells();
		ct.setup_stdcells_mem();

		if (mode_nomux) {
			ct.cell_types.erase(ID($mux));
			ct.cell_types.erase(ID($pmux));
		}

		ct.cell_types.erase(ID($tribuf));
		ct.cell_types.erase(ID($_TBUF_));
		ct.cell_types.erase(ID($anyseq));
		ct.cell_types.erase(ID($anyconst));
		ct.cell_types.erase(ID($allseq));
		ct.cell_types.erase(ID($allconst));

		log("Finding identical cells in module `%s'.\n", module->name);

		int threads = ThreadPool::pool_size(0, module->cells_size());
		int workers = std::max(1, threads);
		OptMergeThreadWorker thread_worker(*module, init_vals, assign_map, ct, workers, mode_share_all, mode_keepdc);

		std::vector<ConcurrentQueue<ComputeCellHashes>> compute_cell_hashes(threads);
		std::vector<ConcurrentQueue<ComputeCellHashesOut>> compute_cell_hashes_out(threads);
		std::vector<ConcurrentQueue<FindDuplicateCells>> find_duplicate_cells(threads);
		std::vector<ConcurrentQueue<FindDuplicateCellsOut>> find_duplicate_cells_out(threads);

		ThreadPool thread_pool(threads, [&](int i) {
			while (std::optional<ComputeCellHashes> c = compute_cell_hashes[i].pop_front()) {
				compute_cell_hashes_out[i].push_back(thread_worker.compute_cell_hashes(*c));
				std::optional<FindDuplicateCells> f = find_duplicate_cells[i].pop_front();
				find_duplicate_cells_out[i].push_back(thread_worker.find_duplicate_cells(i, *f));
			}
		});

		bool did_something = true;
		// A cell may have to go through a lot of collisions if the hash
		// function is performing poorly, but it's a symptom of something bad
		// beyond the user's control.
		while (did_something)
		{
			log("Computing hashes of the cells of `%s'.\n", module->name);
			std::vector<std::vector<std::vector<CellHash>>> bucketed_cell_hashes(workers);
			int cell_index = 0;
			int cells_size = module->cells_size();
			int cells_size_mod_workers = cells_size % workers;
			for (int i = 0; i < workers; ++i) {
				int num_cells = cells_size/workers + ((i < cells_size_mod_workers) ? 1 : 0);
				ComputeCellHashes c = { cell_index, cell_index + num_cells };
				cell_index += num_cells;
				if (threads > 0)
					compute_cell_hashes[i].push_back(c);
				else
					bucketed_cell_hashes[i] = std::move(thread_worker.compute_cell_hashes(c).bucketed_cell_hashes);
			}
			log_assert(cell_index == cells_size);
			if (threads > 0)
				for (int i = 0; i < workers; ++i)
					bucketed_cell_hashes[i] = std::move(compute_cell_hashes_out[i].pop_front()->bucketed_cell_hashes);

			log("Finding duplicate cells in `%s'.\n", module->name);
			std::vector<DuplicateCell> duplicates;
			for (int i = 0; i < workers; ++i) {
				FindDuplicateCells f = { bucketed_cell_hashes };
				if (threads > 0)
					find_duplicate_cells[i].push_back(f);
				else {
					std::vector<DuplicateCell> d = std::move(thread_worker.find_duplicate_cells(i, f).duplicates);
					duplicates.insert(duplicates.end(), d.begin(), d.end());
				}
			}
			if (threads > 0)
				for (int i = 0; i < workers; ++i) {
					std::vector<DuplicateCell> d = std::move(find_duplicate_cells_out[i].pop_front()->duplicates);
					duplicates.insert(duplicates.end(), d.begin(), d.end());
				}
			std::sort(duplicates.begin(), duplicates.end(), [](const DuplicateCell &lhs, const DuplicateCell &rhs) {
				// Sort them by the order in which duplicates would have been detected in a single-threaded
				// run. The cell at which the duplicate would have been detected is the later of the two
				// cells involved.
				return std::max(lhs.remove_cell, lhs.keep_cell) < std::max(rhs.remove_cell, rhs.keep_cell);
			});

			// Convert to cell pointers because removing cells will invalidate the indices.
			std::vector<std::pair<RTLIL::Cell*, RTLIL::Cell*>> cell_ptrs;
			for (DuplicateCell dup : duplicates)
				cell_ptrs.push_back({module->cell_at(dup.remove_cell), module->cell_at(dup.keep_cell)});

			for (auto [remove_cell, keep_cell] : cell_ptrs)
			{
				log_debug("  Cell `%s' is identical to cell `%s'.\n", remove_cell->name, keep_cell->name);
				for (auto &it : remove_cell->connections()) {
					if (remove_cell->output(it.first)) {
						RTLIL::SigSpec keep_sig = keep_cell->getPort(it.first);
						log_debug("    Redirecting output %s: %s = %s\n", it.first,
								log_signal(it.second), log_signal(keep_sig));
						Const init = init_vals(keep_sig);
						init_vals.remove_init(it.second);
						init_vals.remove_init(keep_sig);
						module->connect(RTLIL::SigSig(it.second, keep_sig));
						auto keep_sig_it = keep_sig.begin();
						for (SigBit remove_sig_bit : it.second) {
							assign_map.add(remove_sig_bit, *keep_sig_it);
							++keep_sig_it;
						}
						init_vals.set_init(keep_sig, init);
					}
				}
				log_debug("    Removing %s cell `%s' from module `%s'.\n", remove_cell->type, remove_cell->name, module->name);
				module->remove(remove_cell);
				total_count++;
			}
			did_something = !duplicates.empty();
		}

		for (ConcurrentQueue<ComputeCellHashes> &q : compute_cell_hashes)
			q.close();

		log_suppressed();
	}
};

struct OptMergePass : public Pass {
	OptMergePass() : Pass("opt_merge", "consolidate identical cells") { }
	void help() override
	{
		//   |---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|---v---|
		log("\n");
		log("    opt_merge [options] [selection]\n");
		log("\n");
		log("This pass identifies cells with identical type and input signals. Such cells\n");
		log("are then merged to one cell.\n");
		log("\n");
		log("    -nomux\n");
		log("        Do not merge MUX cells.\n");
		log("\n");
		log("    -share_all\n");
		log("        Operate on all cell types, not just built-in types.\n");
		log("\n");
		log("    -keepdc\n");
		log("        Do not merge flipflops with don't-care bits in their initial value.\n");
		log("\n");
	}
	void execute(std::vector<std::string> args, RTLIL::Design *design) override
	{
		log_header(design, "Executing OPT_MERGE pass (detect identical cells).\n");

		bool mode_nomux = false;
		bool mode_share_all = false;
		bool mode_keepdc = false;

		size_t argidx;
		for (argidx = 1; argidx < args.size(); argidx++) {
			std::string arg = args[argidx];
			if (arg == "-nomux") {
				mode_nomux = true;
				continue;
			}
			if (arg == "-share_all") {
				mode_share_all = true;
				continue;
			}
			if (arg == "-keepdc") {
				mode_keepdc = true;
				continue;
			}
			break;
		}
		extra_args(args, argidx, design);

		int total_count = 0;
		for (auto module : design->selected_modules()) {
			OptMergeWorker worker(module, mode_nomux, mode_share_all, mode_keepdc);
			total_count += worker.total_count;
		}

		if (total_count)
			design->scratchpad_set_bool("opt.did_something", true);
		log("Removed a total of %d cells.\n", total_count);
	}
} OptMergePass;

PRIVATE_NAMESPACE_END
