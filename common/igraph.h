#ifndef CMPLX_COMMON_IGRAPH_H
#define CMPLX_COMMON_IGRAPH_H

#include <igraph.h>
#include <vector>

namespace cmplx {
namespace common {
class IGraph {
 public:
  // nei - the distance (number of steps) within which two vertices will be
  // connected
  static IGraph UndirectedLattice(const std::vector<int>& dimensions, int nei);

  ~IGraph() { igraph_destroy(&graph_); }
  const igraph_t& graph() const { return graph_; }

 private:
  IGraph(const igraph_t& graph_) : graph_(std::move(graph_)) {}
  igraph_t graph_;
};

}  // namespace common
}  // namespace cmplx

#endif  // CMPLX_COMMON_IGRAPH_H
