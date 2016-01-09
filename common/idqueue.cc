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

} // namespace common
} // namespace cmplx

