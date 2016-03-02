#include "./bit_array.h"

#include <algorithm>
#include <cassert>
#include <sstream>

namespace cmplx {
namespace common {
BitArray BitArray::ones(int bits_num) { return ~BitArray(bits_num); }
BitArray BitArray::zeros(int bits_num) { return BitArray(bits_num); }

BitArray::BitArray(int bits_num) : bits_num_(bits_num) {
  bits_in_int_ = 8 * (int)sizeof(uint64_t);
  int vector_size = (bits_num / bits_in_int_) + 1;
  bits_.assign(vector_size, 0);
}

BitArray::BitArray(const BitArray &bit_array) {
  bits_in_int_ = 8 * (int)sizeof(uint64_t);
  bits_num_ = bit_array.bits_num();
  auto data = bit_array.data();
  bits_.assign(data.begin(), data.end());
}

void BitArray::set(int bit_idx, bool value = true) {
  int vector_idx = bit_idx / bits_in_int_;
  bit_idx = bit_idx % bits_in_int_;

  if (value == true) {
    bits_[vector_idx] |= (1LL << bit_idx);
  } else { // value == false
    bits_[vector_idx] &= ~(1LL << bit_idx);
  }
}

bool BitArray::bit(int bit_idx) const {
  int vector_idx = bit_idx / bits_in_int_;
  bit_idx = bit_idx % bits_in_int_;
  return (bits_[vector_idx] >> bit_idx) & 1;
}

int BitArray::bitCount() const {
  int bit_cnt = 0;
  for (uint64_t v : bits_) {
    while (v) {
      v &= (v - 1);
      bit_cnt++;
    }
  }
  return bit_cnt;
}

void BitArray::setWord(int idx, uint64_t value) { bits_[idx] = value; }

std::string BitArray::to_string() const {
  std::ostringstream oss;
  oss << *this;
  return oss.str();
}

/* TODO  Optimize! */
BitArray BitArray::operator~() const {
  int n = bits_num_;
  BitArray ba(n);
  for (int i = 0; i < n; ++i) {
    ba.set(i, ~bit(i));
  }
  return ba;
}

std::ostream &operator<<(std::ostream &os, const BitArray &bit_array) {
  int n = bit_array.bits_num();
  for (int i = 0; i < n; ++i) {
    os << (bit_array.bit(i) ? "1" : "O");
  }
  return os;
}

bool operator==(const BitArray &a, const BitArray &b) {
  if (a.bits_num() != b.bits_num()) {
    return false;
  }

  int data_size = (int)a.data().size();
  for (int i = 0; i < (int)data_size; ++i) {
    if (a.getWord(i) ^ b.getWord(i))
      return false;
  }
  return true;
}

BitArray operator|(const BitArray &a, const BitArray &b) {
  assert(a.bits_num() == b.bits_num());
  BitArray ret(a.bits_num());
  int words = (int)a.data().size();
  for (int i = 0; i < words; ++i) {
    ret.setWord(i, a.getWord(i) | b.getWord(i));
  }
  return ret;
}

BitArray operator&(const BitArray &a, const BitArray &b) {
  assert(a.bits_num() == b.bits_num());
  BitArray ret(a.bits_num());
  int words = (int)a.data().size();
  for (int i = 0; i < words; ++i) {
    ret.setWord(i, a.getWord(i) & b.getWord(i));
  }
  return ret;
}

BitArray operator%(const BitArray &a, const BitArray &b) {
  assert(a.bits_num() == b.bits_num());
  BitArray ret(a.bits_num());
  int words = (int)a.data().size();
  for (int i = 0; i < words - 1; ++i) {
    ret.setWord(i, ~(a.getWord(i) ^ b.getWord(i)));
  }

  int rest = a.bits_num() % (int)sizeof(uint64_t);
  ret.setWord(words - 1, (~(a.getWord(words - 1) ^ b.getWord(words - 1))) &
                             ((1 << rest) - 1));
  return ret;
}

int XNORSimilarity(const BitArray &a, const BitArray &b) {
  BitArray xnor_ba = a % b;
  return xnor_ba.bitCount();
}

std::pair<int, int> JacardSimilarity(const BitArray &a, const BitArray &b) {
  BitArray and_ba = a & b;
  BitArray or_ba = a | b;
  return std::make_pair(and_ba.bitCount(), or_ba.bitCount());
}

} // namespace common
} // namespace cmplx
