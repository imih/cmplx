#include "igraph.h"

#include <cstdio>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "ivector.cc"

namespace cmplx {
namespace common {
IGraph *IGraph::UndirectedLattice(const std::vector<int> dimensions) {
  IVector<int> dim(dimensions);
  igraph_t *graph = (igraph_t *)malloc(sizeof(igraph_t));
  igraph_lattice(graph, dim.vector(), 1, IGRAPH_UNDIRECTED, 0 /* mutual */,
                 0 /* circular */);
  return new IGraph(graph);
}

IGraph *IGraph::BarabasiAlbert(int nodes) {
  igraph_t *graph = (igraph_t *)malloc(sizeof(igraph_t));
  igraph_barabasi_game(graph, nodes, 1, 2, NULL, 0, 1, IGRAPH_UNDIRECTED,
                       IGRAPH_BARABASI_BAG, NULL);
  return new IGraph(graph);
}

IGraph *IGraph::ErdosRenyi(int nodes, double p) {
  igraph_t *graph = (igraph_t *)malloc(sizeof(igraph_t));
  igraph_erdos_renyi_game(graph, IGRAPH_ERDOS_RENYI_GNP, nodes, p,
                          IGRAPH_UNDIRECTED, 0 /* loops*/);
  return new IGraph(graph);
}

/*
bool IGraph::is_connected() const {
  igraph_bool_t bool_t;
  igraph_is_connected(graph_, &bool_t, IGRAPH_WEAK);
  return (bool) bool_t;
}

int IGraph::kCore(int node_id) const {
  if (cores_.get()->empty()) {
    igraph_coreness(graph_, cores_.get()->vector(),
*/
//                    IGRAPH_ALL /* undirected */);
/*
  }
  return (*cores_.get())[node_id];
}

double IGraph::closeness(int node_id) const {
  if (closeness_.get()->empty()) {
    igraph_vs_t *vs_t = (igraph_vs_t *)malloc(sizeof(igraph_vs_t));
    igraph_vs_all(vs_t);
    igraph_closeness(graph_, closeness_.get()->vector(), *vs_t, IGRAPH_ALL,
*/
 //                    NULL, true /* normalize*/);
/*
    delete vs_t;
  }
  return (*closeness_.get())[node_id];
}

double IGraph::betweenness(int node_id) const {
  if (betweenness_.get()->empty()) {
    igraph_vs_t *vs_t = (igraph_vs_t *)malloc(sizeof(igraph_vs_t));
    igraph_vs_all(vs_t);
    igraph_betweenness(graph_, betweenness_.get()->vector(), *vs_t, 0, NULL, 1);
    delete vs_t;
  }
  return (*betweenness_.get())[node_id];
}

void IGraph::writeGML(std::string file_name) {
  FILE* f = fopen(file_name.c_str(), "w");
  igraph_write_graph_gml(graph_, f, NULL, "Iva");
  fclose(f);
}

double IGraph::eigenvector_centrality(int node_id) const {
  if (eigenvector_centrality_.get()->empty()) {
    igraph_arpack_options_t *options =
        (igraph_arpack_options_t *)malloc(sizeof(igraph_arpack_options_t));
    igraph_arpack_options_init(options);
    igraph_eigenvector_centrality(
*/
 //       graph_, eigenvector_centrality_.get()->vector(), NULL, 0 /* directed */,
  //      1 /* scale */, NULL, options);
/*
    delete options;
  }
  return (*eigenvector_centrality_.get())[node_id];
}
*/

IGraph *IGraph::GraphFromGML(const std::string &file_name) {
  igraph_t *graph = (igraph_t *)malloc(sizeof(igraph_t));
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
  return new IGraph(graph);
}

IGraph *IGraph::GraphFromGDF(const std::string &file_name) {
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
  igraph_t *graph = (igraph_t *)malloc(sizeof(igraph_t));
  int err = igraph_create(graph, iedges.vector(), nodes, 0);
  if (err) {
    std::cout << std::string(igraph_strerror(err)) << std::endl;
  }
  return new IGraph(graph);
}

int IGraph::diameter() const {
  igraph_integer_t diam;
  igraph_diameter(graph_, &diam, 0, 0, 0, IGRAPH_UNDIRECTED, 0);
  return (int)diam;
}

const IVector<int> &IGraph::adj_list(int node_id) const {
  if (!adj_list_cache_->count(node_id)) {
    igraph_neighbors(graph_, (*adj_list_cache_)[node_id].vector(), node_id,
                     IGRAPH_OUT);
  igraph_vector_shuffle((*adj_list_cache_)[node_id].vector());
  }
  return adj_list_cache_->at(node_id);
}

}  // namespace common
}  // namespace cmplx
