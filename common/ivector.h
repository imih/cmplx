#ifndef CMPLX_COMMON_IVECTOR_H
#define CMPLX_COMMON_IVECTOR_H

#include <igraph/igraph.h>

#include <vector>

namespace cmplx {
namespace common {
/*
 * Template clas encapsulating igraph_vector_t that hosts igraph_real_t by
 * default.
 * T can be scalar and will be saved as igraph_real_t internally.
 * Purpose: taking care of memory management.
 */
template <class T>
class IVector {
 public:
  /*Initializes an ampty igraph_vector_t.*/
  IVector();
  IVector(const std::vector<T> &v);
  IVector(std::initializer_list<T> il);
  IVector(const IVector<T> &iv);

  ~IVector();

  igraph_vector_t *vector() const { return vector_; }

  T operator[](int idx) const;
  void set(int idx, T val);

  long int size() const { return igraph_vector_size(vector_); }
  bool empty() const { return igraph_vector_empty(vector_); }

 private:
  // By default of type igraph_real_t.
  igraph_vector_t *vector_;
};

template class IVector<int>;

}  // namespace common
}  // namespace cmplx

#endif  // CMPLX_COMMON_IVECTOR_H
