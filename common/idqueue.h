#ifndef CMPLX_COMMON_IDQUEUE_H
#define CMPLX_COMMON_IDQUEUE_H

#include "bit_array.h"

#include <igraph/igraph.h>

#include <cassert>
#include <vector>

namespace cmplx {
namespace common {
/*
 * Class encapsulating igraph_dqueue_t.
 */
class IDqueue {
public:
  IDqueue(long int sz) { assert(!igraph_dqueue_init(&dqueue_, sz)); }

  /* Adds to dqueue all items that are in the bit set.*/
  IDqueue(const BitArray &bit_array);
  ~IDqueue() { igraph_dqueue_destroy(&dqueue_); }

  long int size() { return igraph_dqueue_size(&dqueue_); }
  bool empty() { return igraph_dqueue_empty(&dqueue_); }

  void push(int x) { assert(!igraph_dqueue_push(&dqueue_, x)); }

  int pop() { return (int)igraph_dqueue_pop(&dqueue_); }

private:
  igraph_dqueue_t dqueue_;
};
} // namespace common
} // namespace cmplx

#endif // CMPLX_COMMON_IDQUEUE_H
