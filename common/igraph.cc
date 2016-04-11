#include "igraph.h"

#include <cassert>
#include <cstdio>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>

#include "ivector.h"
#include "ivector.cc"

namespace cmplx {
namespace common {
IGraph IGraph::UndirectedLattice(const std::vector<int> &dimensions) {
  IVector<int> dim(dimensions);
  igraph_t *graph = new igraph_t;
  assert(!igraph_lattice(graph, dim.vector(), 1, IGRAPH_UNDIRECTED,
                         0 /* mutual */, 0 /* circular */));
  return IGraph(graph);
}

IGraph IGraph::GraphFromGML(const std::string &file_name) {
  igraph_t *graph = new igraph_t;
  FILE *f = fopen(file_name.c_str(), "r");
  std::cout << file_name << std::endl;
  if (f == NULL) {
    std::cout << std::strerror(errno) << std::endl;
  }
  int err = igraph_read_graph_gml(graph, f);
  if (err) {
    std::cout << std::string(igraph_strerror(err)) << std::endl;
  }
  fclose(f);
  return IGraph(graph);
}

IGraph IGraph::GraphFromGDF(const std::string &file_name) {
  std::ifstream fs;
  fs.open(file_name);
  if (!fs.is_open()) {
    std::cout << std::strerror(errno) << std::endl;
    exit(1);
  }
  std::string line;
  int nodes;
  std::vector<int> edges;
  while (getline(fs, line)) {
    if (line[0] == 'n' || line[0] == 'e') continue;
    if (line.find(",") == std::string::npos) {
      nodes = stoi(line);
    } else {
      std::istringstream iss(line);
      std::string item;
      while (std::getline(iss, item, ',')) {
        edges.push_back(stoi(item));
      }
    }
  }
  fs.close();

  IVector<int> iedges(edges);
  igraph_t *graph = new igraph_t;
  int err = igraph_create(graph, iedges.vector(), nodes, 0);
  if (err) {
    std::cout << std::string(igraph_strerror(err)) << std::endl;
  }
  return IGraph(graph);
}

int IGraph::diameter() const {
  igraph_integer_t diam;
  igraph_diameter(graph_, &diam, 0, 0, 0, IGRAPH_UNDIRECTED, 0);
  return (int)diam;
}

const IVector<int> &IGraph::adj_list(int node_id) const {
  if (!adj_list_cache_->count(node_id)) {
    assert(!igraph_neighbors(graph_, (*adj_list_cache_)[node_id].vector(),
                             node_id, IGRAPH_OUT));
    igraph_vector_shuffle((*adj_list_cache_)[node_id].vector());
  }
  return adj_list_cache_->at(node_id);
}

}  // namespace common
}  // namespace cmplx
