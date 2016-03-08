#ifndef CMPLX_COMMON_IGRAPH_H
#define CMPLX_COMMON_IGRAPH_H

#include <igraph/igraph.h>
#include <vector>
#include <string>
#include <unordered_map>

#include "ivector.h"

namespace cmplx {
namespace common {

// Represents an undirected graph  that's held in an igraph_t structure.
class IGraph {
public:
  // nei - the distance (number of steps) within which two vertices will be
  // connected
  static IGraph UndirectedLattice(const std::vector<int> &dimensions);
  static IGraph GraphFromGML(const std::string& file_name);
  static IGraph GraphFromGDF(const std::string& file_name);

  IGraph(const IGraph &i_graph) {
    igraph_copy(&graph_, &i_graph.graph());
    adj_list_cache_.clear();
  }

  ~IGraph() { igraph_destroy(&graph_); }

  const igraph_t &graph() const { return graph_; }

  int diameter() const;

  int vertices() const { return (int)igraph_vcount(&graph_); }

  const IVector<int> &adj_list(int node_id) const;

private:
  IGraph(const igraph_t &graph_) : graph_(std::move(graph_)) {}

  igraph_t graph_;
  mutable std::unordered_map<int, IVector<int>> adj_list_cache_;
};

} // namespace common
} // namespace cmplx

#endif // CMPLX_COMMON_IGRAPH_H
