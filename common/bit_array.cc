#include "bit_array.h"

#include <algorithm>

namespace cmplx {
namespace common {
BitArray::BitArray(int bits_num) : bits_num_(bits_num) {
  int vector_size = (bits_num / kFilledBitsNo) + (bits_num % kFilledBitsNo > 0);
  bits_.assign(vector_size, 0);
}

void BitArray::set(int bit_idx, bool value = true) {
  int vector_idx = bit_idx / kFilledBitsNo;
  bit_idx = bit_idx % vector_idx;
  if (value == true) {
    bits_[vector_idx] |= (1LL << bit_idx);
  } else {  // value == false
    bits_[vector_idx] &= ~(1LL << bit_idx);
  }
}

bool BitArray::bit(int bit_idx) const {
  int vector_idx = bit_idx / kFilledBitsNo;
  bit_idx = bit_idx % vector_idx;
  return (bits_[vector_idx] >> bit_idx) & 1;
}

int BitArray::bitCount() const {
  int bit_cnt = 0;
  for (uint64_t v : bits_) {
    // Brian Kernighan bits counting
    while (v != 0) {
      v &= (v - 1);
      bit_cnt++;
    }
  }
  return bit_cnt;
}

int BitArray::XNORSimilarity(const BitArray& bit_array_other) const {
  int bits_num = std::max(bits_num_, bit_array_other.bits_num());
  int bit_cnt = 0;
  for (int i = 0; i < bits_num; ++i) {
    if (bit(i) == bit_array_other.bit(i)) bit_cnt++;
  }
  return bit_cnt;
}

std::pair<int, int> BitArray::JacardSimilarity(
    const BitArray& bit_array_other) const {
  int inter_cnt = 0;
  int union_cnt = 0;
  int bits_num = std::max(bits_num_, bit_array_other.bits_num());
  for (int i = 0; i < bits_num; ++i) {
    bool b1 = bit(i);
    bool b2 = bit_array_other.bit(i);
    if (b1 && b1 == b2) {
      inter_cnt++;
    }
    if (b1 || b2) {
      union_cnt++;
    }
  }
  return std::make_pair(inter_cnt, union_cnt);
}

std::ostream& operator<<(std::ostream& os, const BitArray& bit_array) {
  int n = bit_array.bits_num();
  for (int i = 0; i < n; ++i) {
    os << (bit_array.bit(i) ? "1" : "O");
  }
  return os;
}
}  // namespace common
}  // namespace cmplx

