#include "idqueue.h"

namespace cmplx {
namespace common {
/* Adds to dqueue all items that are in the bit set.*/
IDqueue::IDqueue(const BitArray &bit_array) : IDqueue(bit_array.bits_num()) {
  int sz = bit_array.bits_num();
  for (int i = 0; i < sz; ++i) {
    if (bit_array.bit(i)) {
      push(i);
    }
  }
}

void IDqueue::insertMarked(const BitArray &bit_array) {
  int n = bit_array.bits_num();
  clear();
  for (int i = 0; i < n; ++i) {
    if (bit_array.bit(i))
      push(i);
  }
}

} // namespace common
} // namespace cmplx