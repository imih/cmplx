#include "igraph.h"

#include <cassert>
#include <cstdio>

#include "ivector.h"
#include "ivector.cc"

namespace cmplx {
namespace common {
IGraph IGraph::UndirectedLattice(const std::vector<int>& dimensions) {
  IVector<int> dim(dimensions);
  igraph_t graph;
  assert(!igraph_lattice(&graph, dim.vector(), 1 , IGRAPH_UNDIRECTED,
                         0 /* mutual */, 0 /* circular */));
  return IGraph(std::move(graph));
}

}  // namespace common
}  // namespace cmplx

