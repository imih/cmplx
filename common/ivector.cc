#include "ivector.h"

namespace cmplx {
namespace common {
IVector::IVector(const std::vector<int>& v) {
  assert(igraph_vector_t(&vector_, (int long)v.size()) == 0);
  for (int idx = 0; idx < (int)v.size(); ++idx) {
    VECTOR(vector_)[idx] = v[idx];
  }
}
} // namespace cmplx
} // namespace common

