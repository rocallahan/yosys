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

#include "kernel/yosys.h"
#include "kernel/celltypes.h"
#include "kernel/sigtools.h"
#include "kernel/trtlil.h"

#include <algorithm>

YOSYS_NAMESPACE_BEGIN

TRTLIL::TSigMap TRTLIL::TSigMap::build(const RTLIL::Module &module) {
	TSigMap result;
	int expected_bits = 0;
	for (auto &it : module.connections())
		expected_bits += GetSize(it.first);
	result.bits.reserve(expected_bits*2);
	result.bit_to_net.reserve(expected_bits*2);
	for (auto &it : module.connections()) {
		auto second_it = it.second.begin();
		for (SigBit first_bit : it.first) {
			result.add(first_bit, *second_it);
			++second_it;
		}
	}
	return result;
}

void TRTLIL::TSigMap::add(TSigBit a, TSigBit b) {
	if (a == b)
		return;
	// Put constant in `a` if either of them is a constant.
	if (b.is_const())
		std::swap(a, b);
	int a_index = bits(a);
	int b_index = bits(b);
	bit_to_net.resize(bits.size());
	std::vector<int> *a_net = bit_to_net[a_index];
	std::vector<int> *b_net = bit_to_net[b_index];
	if (a_net == nullptr) {
		if (b_net == nullptr) {
			std::vector<int> *net = new std::vector<int>({a_index, b_index});
			bit_to_net[a_index] = net;
			bit_to_net[b_index] = net;
			return;
		}
		b_net->push_back(a_index);
		if (a.is_const()) {
			std::swap(b_net->front(), b_net->back());
		}
		bit_to_net[a_index] = b_net;
		return;
	}
	if (b_net == nullptr) {
		a_net->push_back(b_index);
		bit_to_net[b_index] = a_net;
		return;
	}
	if (a_net == b_net)
		return;

	// Merge a_net and b_net.
	int a_net_size = GetSize(*a_net);
	int b_net_size = GetSize(*b_net);
	// Put the larger net in `a_net`.
	if (a_net_size < b_net_size) {
		std::swap(a_net, b_net);
		std::swap(a_net_size, b_net_size);
	}
	if (bits[b_net->front()].is_const()) {
		bit_to_net[a_net->front()] = b_net;
		bit_to_net[b_net->front()] = a_net;
		std::swap(a_net->front(), b_net->front());
	}
	a_net->insert(a_net->end(), b_net->begin(), b_net->end());
	for (int bit_index : *b_net)
		bit_to_net[bit_index] = a_net;
	delete b_net;
}

TRTLIL::TSigMap::~TSigMap() {
	for (int i = 0; i < GetSize(bit_to_net); ++i) {
		std::vector<int> *net = bit_to_net[i];
		if (net == nullptr)
			continue;
		for (int bit_index : *net) {
			bit_to_net[bit_index] = nullptr;
		}
		delete net;
	}
}

TRTLIL::TSigSpec TRTLIL::TSigSpec::build(const RTLIL::Const &c) {
	std::vector<TSigBit> bits;
	bits.reserve(c.size());
	for (RTLIL::State s : c) {
		bits.emplace_back(s);
	}
	return std::move(bits);
}

TRTLIL::TSigSpec TRTLIL::TSigSpec::build(const RTLIL::SigSpec &sigspec) {
	std::vector<TSigBit> bits;
	bits.reserve(sigspec.size());
	for (RTLIL::SigBit bit : sigspec) {
		bits.emplace_back(bit);
	}
	return std::move(bits);
}

TRTLIL::TSigSpec TRTLIL::TSigSpec::extract(int offset, int length) const {
	TSigSpec result;
	result.bits.insert(result.bits.end(), bits.begin() + offset, bits.begin() + offset + length);
	return result;
}

void TRTLIL::TSigSpec::sort() {
	std::sort(bits.begin(), bits.end());
}

void TRTLIL::TSigSpec::sort_and_unify() {
	std::sort(bits.begin(), bits.end());
	auto last = std::unique(bits.begin(), bits.end());
	bits.erase(last, bits.end());
}

YOSYS_NAMESPACE_END
