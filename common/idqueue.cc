#include "idqueue.h"

namespace cmplx {
namespace common {
IDqueue::IDqueue(long int sz) {
  dqueue_ = (igraph_dqueue_t*) malloc(sizeof(igraph_dqueue_t));
  igraph_dqueue_init(dqueue_, sz);
  clear();
}

/* Adds to dqueue all items that are in the bit set.*/
IDqueue::IDqueue(const BitArray &bit_array) : IDqueue(bit_array.bits_num()) {
  int sz = bit_array.bits_num();
  for (int i = 0; i < sz; ++i) {
    if (bit_array.bit(i)) {
      push(i);
    }
  }
}

IDqueue::~IDqueue() {
  igraph_dqueue_destroy(dqueue_);
  if(dqueue_) igraph_free(dqueue_);
}

void IDqueue::insertMarked(const BitArray &bit_array) {
  int n = bit_array.bits_num();
  clear();
  for (int i = 0; i < n; ++i) {
    if (bit_array.bit(i)) push(i);
  }
}

}  // namespace common
}  // namespace cmplx
