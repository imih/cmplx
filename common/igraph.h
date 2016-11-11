#ifndef CMPLX_COMMON_IGRAPH_H
#define CMPLX_COMMON_IGRAPH_H

#include <igraph/igraph.h>
#include <vector>
#include <string>
#include <unordered_map>

#include "ivector.h"

#include <memory>

namespace cmplx {
namespace common {

/**
 * Represents an undirected graph  that's held in an igraph_t structure.
 */
class IGraph {
 public:
  static IGraph* UndirectedLattice(const std::vector<int> dimensions);
  static IGraph* GraphFromGML(const std::string& file_name);
  static IGraph* GraphFromGDF(const std::string& file_name);
  static IGraph* BarabasiAlbert(int nodes);
  static IGraph* ErdosRenyi(int nodes, double p);

  IGraph(const IGraph& i_graph);
  ~IGraph();

  const igraph_t& graph() const { return *graph_; }

  int diameter() const;
  bool is_connected() const;

  int vertices() const { return (int)igraph_vcount(graph_); }

  const IVector<int>& adj_list(int node_id) const;

  int deg(int node_id) const { return adj_list(node_id).size(); }
  int kCore(int node_id) const;
  double closeness(int node_id) const;
  double betweenness(int node_id) const;
  double eigenvector_centrality(int node_id) const;

  void writeGML(std::string file_name);

 private:
  IGraph(igraph_t* graph_, bool fill_attributes = false);

  igraph_t* graph_;
  mutable std::unordered_map<int, IVector<int> >* adj_list_cache_;

  mutable std::unique_ptr<IVector<int> > cores_;
  mutable std::unique_ptr<IVector<double> > closeness_;
  mutable std::unique_ptr<IVector<double> > betweenness_;
  mutable std::unique_ptr<IVector<double> > eigenvector_centrality_;
};

}  // namespace common
}  // namespace cmplx

#endif  // CMPLX_COMMON_IGRAPH_H
