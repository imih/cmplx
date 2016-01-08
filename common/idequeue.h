#ifndef CMPLX_COMMON_IDQUEUE_H
#define CMPLX_COMMON_IDQUEUE_H

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
  IDqueue(long int size) { assert(!igraph_dqueue_init(&dqueue_, size)); }

  ~IDqueue() { igraph_dqueue_destroy(&dqueue_); }

  long int size() { return igraph_dqueue_size(&dqueue_); }
  bool empty() { return igraph_dqueue_empty(&dqueue_); }

 private:
  igraph_dqueue_t dqueue_;
};
}  // namespace common
}  // namespace cmplx

#endif  // CMPLX_COMMON_IDQUEUE_H
