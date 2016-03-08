#include "direct_mc_params.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::Realization;

namespace {
std::vector<std::string> split(std::string s) {
  std::stringstream ss(s);
  std::vector<std::string> elems;
  std::string item;
  while (std::getline(ss, item, ' ')) {
    elems.push_back(item);
  }
  return elems;
}

const int SOURCE_LINE = 0;
const int P_LINE = 1;
const int Q_LINE = 2;
const int T_LINE = 3;
const int NODES_LINE = 5;
const int SIMUL_LINE = 0;
std::string BENCHMARK_PATH =
    "/home/iva/dipl/nino-sup/Supplementary_data_code/Data/BenchmarkData";
}

namespace cmplx {
DirectMCParams DirectMCParams::SupFig2Params() {
  int lattice_size1 = 5;
  int lattice_size2 = 4;
  IGraph graph = IGraph::UndirectedLattice({lattice_size1, lattice_size2});

  int vertices = graph.vertices();
  BitArray r = BitArray::zeros(vertices);
  std::vector<int> infected = {2, 6, 7, 8, 9, 12};
  for (int inf_id : infected) {
    r.set(inf_id, true);
  }
  double p = 0.2;
  double q = 0.3;
  int maxT = 5;
  Realization realization = Realization(p, q, maxT, r);

  int simulations = 1000000000;
  return DirectMCParams(graph, realization, simulations);
}

// TODO determine number of simulations yourself!
DirectMCParams DirectMCParams::BenchmarkParams(int realization_no) {
  IGraph graph = IGraph::GraphFromGDF(BENCHMARK_PATH + "/network/lattice_gephi.GDF");
  double p = 0, q = 0;
  int T = 0;
  BitArray r = BitArray::zeros(graph.vertices());

  std::ifstream f_real;
  f_real.open(BENCHMARK_PATH + "/realizations/realization_" +
              std::to_string(realization_no) + ".txt");
  if (!f_real.is_open()) {
    std::cout << std::strerror(errno) << std::endl;
    exit(1);
  }
  std::string line;
  int line_no = 0;
  while (getline(f_real, line)) {
    auto items = split(line);
    switch (line_no) {
    case P_LINE:
      p = stod(items[1]);
      break;
    case Q_LINE:
      q = stod(items[1]);
      break;
    case T_LINE:
      T = stoi(items[1]);
      break;
    }

    if (line_no >= NODES_LINE) {
      if (stoi(line)) {
        r.set(line_no - NODES_LINE, true);
      }
    }
    line_no++;
  }

  f_real.close();

  std::ifstream f_sol;
  f_sol.open(BENCHMARK_PATH + "/solutions/inverse_solution_" +
             std::to_string(realization_no) + ".txt");
  if (!f_sol.is_open()) {
    std::cout << std::strerror(errno) << std::endl;
    exit(1);
  }
  int simulations = 0;
  getline(f_sol, line);
  f_sol.close();
  auto items = split(line);
  simulations = stoi(items[3]);

  Realization realization(p, q, T, r);
  return DirectMCParams(graph, realization, simulations);
}

} // namespace cmplx

