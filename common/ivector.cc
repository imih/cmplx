#include "ivector.h"

#include <igraph/igraph.h>
#include <cassert>

namespace cmplx {
namespace common {
template <class T>
IVector<T>::IVector(const std::vector<T>& v) {
  assert(igraph_vector_init(&vector_, (int long)v.size()) == 0);
  for (int idx = 0; idx < (int)v.size(); ++idx) {
    VECTOR(vector_)[idx] = v[idx];
  }
}

template <class T>
IVector<T>::IVector(std::initializer_list<T> il) {
  assert(igraph_vector_init(&vector_, (int long)il.size()) == 0);
  int idx = 0;
  for (T v : il) {
    VECTOR(vector_)[idx] = v;
    idx++;
  }
}

template <class T>
T& IVector<T>::operator[](int idx) {
  assert(idx >= 0);
  assert(idx < (int)size());
  return (T&) VECTOR(vector_)[idx];
}
}  // namespace cmplx
}  // namespace common

