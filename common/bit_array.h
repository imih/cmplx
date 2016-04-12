#ifndef CMPLX_COMMON_BIT_ARRAY_H
#define CMPLX_COMMON_BIT_ARRAY_H

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace cmplx {
namespace common {
class BitArray {
public:
  static BitArray ones(int bits_num);
  static BitArray zeros(int bits_num);

  BitArray(int bits_num);
  BitArray(const BitArray &bit_array);
  ~BitArray() = default;

  void set(int bit_idx, bool value);

  bool bit(int bit_idx) const;
  int bits_num() const { return bits_num_; }
  int bitCount() const;

  const std::vector<uint64_t> &data() const { return bits_; }
  void setWord(int idx, uint64_t value);
  uint64_t getWord(int idx) const { return bits_[idx]; }

  std::string to_string() const;

  BitArray operator~() const;

  std::vector<int> positions() const;

private:
  int bits_num_;
  int bits_in_int_;

  std::vector<uint64_t> bits_;
};

std::ostream &operator<<(std::ostream &oc, const BitArray &bit_array);

bool operator==(const BitArray &a, const BitArray &b);

// OR
BitArray operator|(const BitArray &a, const BitArray &b);
// AND
BitArray operator&(const BitArray &a, const BitArray &b);
// XNOR
BitArray operator%(const BitArray &a, const BitArray &b);

int XNORSimilarity(const BitArray &a, const BitArray &b);

double JaccardSimilarity(const BitArray &a, const BitArray &b);

} // namespace common
} // namespace cmplx

#endif // CMPLX_COMMON_BIT_ARRAY_H
