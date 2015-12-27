#ifndef CMPLX_COMMON_IVECTOR_H
#define CMPLX_COMMON_IVECTOR_H

#include <igraph.h>

#include <vector>

namespace cmplx {
namespace common {
class IVector {
 public:
  IVector(const std::vector<int>& v);

  ~IVector() { igraph_vector_destroy(&vector_); }

  const igraph_vector_t& vector() const { return vector_; }

 private:
  igraph_vector_t vector_;
};
}
}

#endif  // CMPLX_COMMON_IVECTOR_H
