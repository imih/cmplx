#ifndef CMPLX_COMMON_BIT_ARRAY_H
#define CMPLX_COMMON_BIT_ARRAY_H

#include <iostream>
#include <vector>
#include <cstdint>

namespace {
const int kFilledBitsNo = 60;
}  // namespace

namespace cmplx {
namespace common {
class BitArray {
 public:
  BitArray(int bits_num);

  void set(int bit_idx, bool value);

  bool bit(int bit_idx) const;
  int bits_num() const { return bits_num_; }

  int bitCount() const;
  // XNOR_bitCount
  int XNORSimilarity(const BitArray& bit_array_other) const;

  // Jacard_similarity is first_num / second_num
  std::pair<int, int> JacardSimilarity(const BitArray& bit_array_other) const;

 private:
  int bits_num_;
  std::vector<uint64_t> bits_;
};

std::ostream& operator<<(std::ostream& oc, const BitArray& bit_array);

}  // namespace common
}  // namespace cmplx

#endif  // CMPLX_COMMON_BIT_ARRAY_H
