#include "ivector.h"

#include <iostream>
#include <igraph/igraph.h>
#include <cassert>

namespace cmplx {
namespace common {
template <class T>
IVector<T>::IVector(const std::vector<T>& v) {
  assert(igraph_vector_init(&vector_, (int long)v.size()) == 0);
  int idx = 0;
  for (T val : v) {
    VECTOR(vector_)[idx] = val;
    idx++;
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
T IVector<T>::operator[](int idx) {
  assert(idx >= 0);
  assert(idx < (int)size());
  return (T)VECTOR(vector_)[idx];
}

template <class T>
void IVector<T>::set(int idx, T val) {
  assert(idx >= 0);
  assert(idx < (int)size());
  VECTOR(vector_)[idx] = val;
}

}  // namespace cmplx
}  // namespace common

