/* -*- c++ -*-
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

#ifndef TRTLIL_H
#define TRTLIL_H

#include "kernel/rtlil.h"
#include "kernel/celltypes.h"
#include "kernel/ffinit.h"

YOSYS_NAMESPACE_BEGIN

/**
 * The R... classes are read-only wrappers around RTLIL objects. They
 * are safe for multithreaded access to RTLIL (as long as no other access to
 * RTLIL occurs concurrently). In particular they are careful not to bump
 * refcounts or trigger representation changes in SigSpec or SigMap.
 *
 * The T... classes are thread-compatible data structures. The const
 * methods are safe for multithreaded access, and the non-const methods
 * are not.
 *
 * Public methods generally don't mention "unsafe" RTLIL types.
 * Safe RTLIL types are limited to:
 * -- State
 * -- Const
 * -- StaticIdString
 * Constructors and constructor-like functions are allowed to take unsafe RTLIL types
 * as parameters.
 */
namespace TRTLIL
{
	class RIdString {
	public:
		RIdString(const RTLIL::IdString &str) : str(&str) {}
		RIdString(const RTLIL::StaticIdString &str) : str(&str.operator const Yosys::RTLIL::IdString &()) {}
		RIdString(const RIdString &) = default;
		template<typename... Args>
		bool in(const Args &... args) const {
			return (... || in(args));
		}
		bool in(const RTLIL::StaticIdString &rhs) const { return *this == rhs; }
		bool operator==(const RTLIL::StaticIdString &rhs) const { return *str == rhs; }
		bool operator!=(const RTLIL::StaticIdString &rhs) const { return *str != rhs; }
		bool operator==(const RIdString &rhs) const { return *str == *rhs.str; }
		bool operator!=(const RIdString &rhs) const { return *str != *rhs.str; }

		[[nodiscard]] Hasher hash_into(Hasher h) const { h.eat(*str); return h; }
	private:
		friend class RCellTypes;
		friend class RCell;

		const RTLIL::IdString *str;
	};

	class TSigBit {
	public:
		TSigBit(const RTLIL::SigBit &bit) : bit(bit) {}
		TSigBit(const TSigBit &) = default;
		TSigBit &operator=(const TSigBit &) = default;
		bool operator==(const TSigBit &rhs) const { return bit == rhs.bit; }
		bool operator!=(const TSigBit &rhs) const { return bit != rhs.bit; }
		bool operator<(const TSigBit &rhs) const { return bit < rhs.bit; }

		bool is_const() const { return bit.wire == nullptr; }

		[[nodiscard]] Hasher hash_into(Hasher h) const { h.eat(bit); return h; }
	private:
		friend class RFfInitVals;
		friend class TSigMap;
		RTLIL::SigBit bit;
	};

	class TSigSpec {
	public:
		TSigSpec() {}
		TSigSpec(TSigSpec &&) = default;
		TSigSpec(const TSigSpec &) = default;
		TSigSpec(std::vector<TSigBit> bits) : bits(std::move(bits)) {}
		TSigSpec &operator=(TSigSpec &&) = default;
		TSigSpec &operator=(const TSigSpec &) = default;

		static TSigSpec build(const RTLIL::SigSpec &sigspec);
		static TSigSpec build(const RTLIL::Const &c);

		TSigSpec extract(int offset, int length) const;
		int size() const { return bits.size(); }
		bool operator==(const TSigSpec &rhs) const { return bits == rhs.bits; }
		bool operator!=(const TSigSpec &rhs) const { return bits != rhs.bits; }
		bool operator<(const TSigSpec &rhs) const { return bits < rhs.bits; }

		void sort();
		void sort_and_unify();
		void append(TSigBit bit) { bits.push_back(bit); }
		void append(const TSigSpec &sigspec) { bits.insert(bits.end(), sigspec.bits.begin(), sigspec.bits.end()); }

		[[nodiscard]] Hasher hash_into(Hasher h) const { h.eat(bits); return h; }

		using const_iterator = std::vector<TSigBit>::const_iterator;
		const_iterator begin() const { return bits.begin(); }
		const_iterator end() const { return bits.end(); }
	private:
		std::vector<TSigBit> bits;
	};

	// For thread-safety we can't use union-find to track nets (i.e. equivalence classes
	// of connected signals), because lookup must not modify the data in any way.
	// Instead, we keep a map from SigBit to equivalence class, and for each equivalence
	// class, store the complete list of SigBits mapping to that class.
	// When we merge two equivalence classes, we eagerly move all SigBits from the smaller
	// class into the larger class. This guarantees that a series of N add() operations takes
	// at most O(N log N) time.
	class TSigMap {
	public:
		TSigMap() = default;
		TSigMap(TSigMap &&) = default;
		TSigMap &operator=(TSigMap &&rhs) = default;
		~TSigMap();

		// Call this only on the main thread.
		// Do we need to parallelize this?
		static TSigMap build(const RTLIL::Module &module);

		void add(TSigBit a, TSigBit b);

		TSigBit map(TSigBit bit) const {
			int bit_index = bits.lookup(bit);
			if (bit_index < 0)
				return bit;
			return bits[bit_to_net[bit_index]->front()];
		}
	private:
		friend class RCell;
		friend struct GenericFfInitVals<TSigMap>;

		TSigSpec map(const RTLIL::SigSpec &sigspec) const {
			std::vector<TSigBit> bits;
			bits.reserve(sigspec.size());
			for (SigBit bit : sigspec) {
				bits.push_back((*this)(bit));
			}
			return TSigSpec(std::move(bits));
		}
		SigBit operator()(RTLIL::SigBit bit) const {
			return map(TSigBit(bit)).bit;
		}
		SigSpec operator()(const RTLIL::SigSpec &sigspec) const {
			std::vector<SigBit> bits;
			bits.reserve(sigspec.size());
			for (SigBit bit : sigspec) {
				bits.push_back((*this)(bit));
			}
			return SigSpec(std::move(bits));
		}

		// A bit is only present here if it belongs to a net of size > 1.
		idict<TSigBit> bits;
		// Every net has an allocated std::vector<int> containing the indices of
		// its SigBits. Only nets of size > 1 are stored, so this vector must have
		// at least 2 elements.
		// For every SigBit in the net, there is an entry in this dict with a pointer to
		// that shared std::vector.
		// If one of the SigBits is a constant, then the first element of the vector
		// will be a constant.
		std::vector<std::vector<int>*> bit_to_net;
	};

	class RFfInitVals {
	public:
		RFfInitVals(const GenericFfInitVals<TSigMap> &init_vals) : init_vals(init_vals) {}
		RTLIL::State operator()(TSigBit bit) const { return init_vals(bit.bit); }
		RTLIL::Const operator()(const TSigSpec &sigspec) const {
			RTLIL::Const::Builder builder(sigspec.size());
			for (TSigBit bit : sigspec) {
				builder.push_back((*this)(bit));
			}
			return builder.build();
		}
	private:
		friend class RCell;
		const GenericFfInitVals<TSigMap> &init_vals;
	};

	class RCell {
	public:
		struct mapped_port_iterator {
		public:
			using iterator_category = std::forward_iterator_tag;
			using value_type = std::pair<RIdString, TSigSpec>;
			using difference_type = ptrdiff_t;
			using pointer = value_type*;
			using reference = value_type&;

			value_type operator*() const {
				return {it->first, sigmap.map(it->second)};
			}
			mapped_port_iterator &operator++() { ++it; return *this; }
			bool operator==(const mapped_port_iterator &rhs) const { return it == rhs.it; }
			bool operator!=(const mapped_port_iterator &rhs) const { return it != rhs.it; }
		private:
			friend class RCell;
			mapped_port_iterator(const TSigMap &sigmap, dict<RTLIL::IdString, RTLIL::SigSpec>::const_iterator it)
					: sigmap(sigmap), it(it) {}
			const TSigMap &sigmap;
			dict<RTLIL::IdString, RTLIL::SigSpec>::const_iterator it;
		};

		struct MappedPorts {
		public:
			mapped_port_iterator begin() const { return {sigmap, ports.begin()}; }
			mapped_port_iterator end() const { return {sigmap, ports.end()}; }
			int size() const { return ports.size(); }
		private:
			friend class RCell;
			MappedPorts(const TSigMap &sigmap, const dict<RTLIL::IdString, RTLIL::SigSpec> &ports)
					: sigmap(sigmap), ports(ports) {}
			const TSigMap &sigmap;
			const dict<RTLIL::IdString, RTLIL::SigSpec> &ports;
		};

		struct named_const_iterator {
		public:
			using iterator_category = std::forward_iterator_tag;
			using value_type = std::pair<RIdString, const RTLIL::Const&>;
			using difference_type = ptrdiff_t;
			using pointer = value_type*;
			using reference = value_type&;

			value_type operator*() const {
				const auto &[id, c] = *it;
				return {id, c};
			}
			named_const_iterator &operator++() { ++it; return *this; }
			bool operator==(const named_const_iterator &rhs) const { return it == rhs.it; }
			bool operator!=(const named_const_iterator &rhs) const { return it != rhs.it; }
		private:
			friend class RCell;
			friend class NamedConsts;
			named_const_iterator(dict<RTLIL::IdString, RTLIL::Const>::const_iterator it)
					: it(it) {}
			dict<RTLIL::IdString, RTLIL::Const>::const_iterator it;
		};

		struct NamedConsts {
		public:
			named_const_iterator begin() const { return {consts.begin()}; }
			named_const_iterator end() const { return {consts.end()}; }
			bool operator==(const NamedConsts &rhs) const { return consts == rhs.consts; }
			bool operator!=(const NamedConsts &rhs) const { return consts != rhs.consts; }
			int size() const { return consts.size(); }
		private:
			friend class RCell;
			NamedConsts(const dict<RTLIL::IdString, RTLIL::Const> &consts)
					: consts(consts) {}
			const dict<RTLIL::IdString, RTLIL::Const> &consts;
		};

		RCell(const RTLIL::Cell &cell) : cell(cell) {}

		RIdString get_type() const { return cell.type; }
		bool is_builtin_ff() const { return cell.is_builtin_ff(); }
		bool known() const { return cell.known(); }

		bool has_port(RIdString name) const { return cell.hasPort(*name.str); }
		TSigSpec get_port_mapped(RIdString name, const TSigMap &map) const {
			return map.map(cell.getPort(*name.str));
		}
		RTLIL::Const get_port_init_vals(RIdString name, const RFfInitVals &init_vals) const {
			return init_vals.init_vals(cell.getPort(*name.str));
		}
		MappedPorts mapped_ports(const TSigMap &sigmap) const { return {sigmap, cell.connections()}; }
		bool output(RIdString port_name) const { return cell.output(*port_name.str); }

		NamedConsts parameters() const { return {cell.parameters}; }
		NamedConsts attributes() const { return {cell.attributes}; }
		bool has_keep_attr() const { return cell.has_keep_attr(); }

		bool same_cell(const RCell &rhs) const { return &cell == &rhs.cell; }

	private:
		friend class RCellTypes;
		friend class RModule;

		const RTLIL::Cell &cell;
	};

	class RCellTypes {
	public:
		RCellTypes(const CellTypes &cell_types) : cell_types(cell_types) {}
		bool cell_known(RIdString cell_type) const { return cell_types.cell_known(*cell_type.str); }
	private:
		const CellTypes &cell_types;
	};

	class RModule {
	public:
		RModule(const RTLIL::Module &module) : module(module) {}
		int cells_size() const { return module.cells_size(); }
		RCell cell_at(int index) const { return *module.cell_at(index); }
		bool selected(RCell cell) const {
			return module.design->selected(&module, &cell.cell);
		}
	private:
		const RTLIL::Module &module;
	};
}

YOSYS_NAMESPACE_END

#endif
