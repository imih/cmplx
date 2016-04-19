#include "ivector.h"

#include <iostream>
#include <igraph/igraph.h>
#include <cassert>

namespace cmplx {
namespace common {
template <class T>
IVector<T>::IVector() {
  vector_ = (igraph_vector_t*) malloc(sizeof(igraph_vector_t));
  !igraph_vector_init(vector_, 0);
}

template <class T>
IVector<T>::IVector(const std::vector<T> &v) {
  vector_ = (igraph_vector_t*) malloc(sizeof(igraph_vector_t));
  igraph_vector_init(vector_, (int long)v.size());
  //printf("%ld\n", size());
  for(int i = 0; i < (int) v.size(); ++i) {
    //igraph_vector_set(vector_, i, v[i]);
    VECTOR(*vector_)[i] = v[i];
  }
}

template <class T>
IVector<T>::IVector(std::initializer_list<T> il) {
  vector_ = (igraph_vector_t*) malloc(sizeof(igraph_vector_t));
  igraph_vector_init(vector_, (int long)il.size());
  int idx = 0;
  for (T v : il) {
    VECTOR(*vector_)[idx] = v;
    idx++;
  }
}

template <class T>
T IVector<T>::operator[](int idx) const {
  assert(idx >= 0);
  assert(idx < (int)size());
  return (T)VECTOR(*vector_)[idx];
}

template <class T>
void IVector<T>::set(int idx, T val) {
  assert(idx >= 0);
  assert(idx < (int)size());
  VECTOR(*vector_)[idx] = val;
}

}  // namespace cmplx
}  // namespace common
