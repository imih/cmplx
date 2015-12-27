#include "igraph.h"
#include "ivector.h"

namespace cmplx {
namespace common {
IGraph IGraph::UndirectedLattice(const std::vector<int>& dimensions, int nei) {
  assert(!igraph_lattice(&graph_, IVector(dimensions), nei, IGRAPH_UNDIRECTED,
                         FALSE /* mutual */, FALSE /* circular */));
}

}  // namespace common
}  // namespace cmplx

