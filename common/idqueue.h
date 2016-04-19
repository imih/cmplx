#ifndef CMPLX_COMMON_IDQUEUE_H
#define CMPLX_COMMON_IDQUEUE_H

#include "bit_array.h"

#include <igraph/igraph.h>

#include <vector>

namespace cmplx {
namespace common {
/*
 * Class encapsulating igraph_dqueue_t.
 */
class IDqueue {
 public:
  IDqueue(long int sz);
  /* Adds to dqueue all items that are in the bit set.*/
  IDqueue(const BitArray &bit_array);

  ~IDqueue();

  long int size() { return igraph_dqueue_size(dqueue_); }
  bool empty() { return igraph_dqueue_empty(dqueue_); }

  void push(int x) { igraph_dqueue_push(dqueue_, x); }

  void clear() { if(size()) igraph_dqueue_clear(dqueue_); }

  int pop() { return (int)igraph_dqueue_pop(dqueue_); }

  void insertMarked(const BitArray &bit_array);

 private:
  IDqueue(const IDqueue&);

  igraph_dqueue_t *dqueue_;
};
}  // namespace common
}  // namespace cmplx

#endif  // CMPLX_COMMON_IDQUEUE_H
