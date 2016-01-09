#ifndef CMPLX_COMMON_IGRAPH_H
#define CMPLX_COMMON_IGRAPH_H

#include <igraph/igraph.h>
#include <vector>
#include <unordered_map>

#include "ivector.h"

namespace cmplx {
namespace common {
 
// Represents an undirected graph  that's held in an igraph_t structure.
class IGraph {
 public:
  // nei - the distance (number of steps) within which two vertices will be
  // connected
  static IGraph UndirectedLattice(const std::vector<int>& dimensions);

  IGraph(const IGraph& i_graph) {
    igraph_copy(&graph_, &i_graph.graph());
  }

  ~IGraph() { igraph_destroy(&graph_); }

  const igraph_t& graph() const { return graph_; }

  int diameter() {
    igraph_integer_t diam;
    igraph_diameter(&graph_, &diam, 0, 0, 0, IGRAPH_UNDIRECTED, 0);
    return (int) diam;
  }

  int vertices() {
    return (int)  igraph_vcount(&graph_);
  }

  const IVector<int> &adj_list(int node_id) {
    if (!adj_list_cache_.count(node_id)) {
      IVector<int> v;
      assert(!igraph_incident(&graph_, v.vector(), node_id, IGRAPH_OUT));
      adj_list_cache_[node_id] = v;
    }
    return adj_list_cache_[node_id];
  }

 private:
  IGraph(const igraph_t& graph_) : graph_(std::move(graph_)) {}
  igraph_t graph_;
  std::unordered_map<int, IVector<int> > adj_list_cache_;
};

}  // namespace common
}  // namespace cmplx

#endif  // CMPLX_COMMON_IGRAPH_H
