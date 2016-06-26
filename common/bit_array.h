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

  BitArray() : BitArray(0) {}
  BitArray(int bits_num);
  BitArray(const BitArray &bit_array);
  ~BitArray() = default;

  void set(int bit_idx, bool value);

  bool bit(int bit_idx) const;

  /** Returns number of bits represented by the bitset.
   */
  int bits_num() const { return bits_num_; }

  /** Returns number of bits set to 1.
   */
  int bitCount() const;

  /** Returns array of indices corresponding to bits set to 1.
   */
  std::vector<int> positions() const;

  std::string to_string() const;

  friend std::ostream &operator<<(std::ostream &oc, const BitArray &bit_array);

  friend bool operator==(const BitArray &a, const BitArray &b);

  // OR
  friend BitArray operator|(const BitArray &a, const BitArray &b);
  // AND
  friend BitArray operator&(const BitArray &a, const BitArray &b);
  // XNOR
  friend BitArray operator%(const BitArray &a, const BitArray &b);

  friend int XNORSimilarity(const BitArray &a, const BitArray &b);

  friend double JaccardSimilarity(const BitArray &a, const BitArray &b);

 private:
  void setWord(int idx, uint64_t value);
  uint64_t getWord(int idx) const { return bits_[idx]; }

  BitArray operator~() const;

  int bits_num_;
  int bits_in_int_;

  std::vector<uint64_t> bits_;
};

}  // namespace common
}  // namespace cmplx

#endif  // CMPLX_COMMON_BIT_ARRAY_H
