#include "igraph.h"

#include <cassert>
#include <cstdio>

#include "ivector.h"
#include "ivector.cc"

namespace cmplx {
namespace common {
IGraph IGraph::UndirectedLattice(const std::vector<int> &dimensions) {
  IVector<int> dim(dimensions);
  igraph_t graph;
  assert(!igraph_lattice(&graph, dim.vector(), 1, IGRAPH_UNDIRECTED,
                         0 /* mutual */, 0 /* circular */));
  return IGraph(std::move(graph));
}

int IGraph::diameter() {
  igraph_integer_t diam;
  igraph_diameter(&graph_, &diam, 0, 0, 0, IGRAPH_UNDIRECTED, 0);
  return (int)diam;
}

const IVector<int> &IGraph::adj_list(int node_id) {
  if (!adj_list_cache_.count(node_id)) {
    IVector<int> v;
    assert(!igraph_incident(&graph_, v.vector(), node_id, IGRAPH_OUT));
    adj_list_cache_[node_id] = v;
  }
  return adj_list_cache_[node_id];
}

} // namespace common
} // namespace cmplx

