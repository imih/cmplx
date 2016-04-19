#include "source_detection_params.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <random>

#include "simul/simulator.h"

using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::Realization;
using cmplx::common::SirParams;

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
    "/home/iva/dipl/Supplementary_data_code/Data/Benchmark data";
}

namespace cmplx {
std::unique_ptr<SourceDetectionParams> SourceDetectionParams::SupFig2Params() {
  int lattice_size1 = 5;
  int lattice_size2 = 4;
  std::vector<int> size = std::vector<int>({lattice_size1, lattice_size2});
  IGraph* graph = IGraph::UndirectedLattice(size);

  int vertices = graph->vertices();
  BitArray r = BitArray::zeros(vertices);
  std::vector<int> infected = {2, 6, 7, 8, 9, 12};
  for (int inf_id : infected) {
    r.set(inf_id, true);
  }
  double p = 0.2;
  double q = 0.3;
  int maxT = 5;
  Realization realization =
      Realization(p, q, maxT, r, BitArray::zeros(vertices));

  int simulations = 1000000000; /* supposed to be 10e9 */
  return std::unique_ptr<SourceDetectionParams>(
      new SourceDetectionParams(graph, realization, simulations));
}

std::unique_ptr<SourceDetectionParams> SourceDetectionParams::LatticeCenter() {
  int lattice_size = 3;
  IGraph* graph = IGraph::UndirectedLattice({lattice_size, lattice_size});
  BitArray r = BitArray::ones(graph->vertices());
  double p = 0.2;
  double q = 0;
  Realization realization =
      Realization(p, q, 2, r, BitArray::zeros(graph->vertices()));
  return std::unique_ptr<SourceDetectionParams>(
      new SourceDetectionParams(graph, realization, 1000000));
}

namespace {
int chooseSource(int n) {
  std::uniform_int_distribution<int> d(0, n - 1);
  static thread_local std::mt19937 gen;
  return d(gen);
}
}

std::unique_ptr<SourceDetectionParams> SourceDetectionParams::ParamsFromGrid(
    double p, double q, int n) {
  IGraph* graph = IGraph::UndirectedLattice({n, n});
  int TMax = 5;
  int source_v = n * n / 2;

  cmplx::simul::Simulator simulator(graph);
  while (true) {
    BitArray r = BitArray::zeros(graph->vertices());
    r.set(source_v, true);
    BitArray s = BitArray::ones(graph->vertices());
    s.set(source_v, false);
    SirParams sir_params(p, q, TMax, r, s);
    simulator.NaiveSIR(sir_params);
    Realization real(p, q, TMax, sir_params.infected(), sir_params.recovered());
    if (real.realization().bitCount() > 1)
      return std::unique_ptr<SourceDetectionParams>(
          new SourceDetectionParams(graph, real, 1000000));
  }
}

// TODO determine number of simulations yourself!
std::unique_ptr<SourceDetectionParams> SourceDetectionParams::BenchmarkParams(
    int realization_no) {
  IGraph* graph =
      IGraph::GraphFromGML(BENCHMARK_PATH + "/network/lattice5.gml");
  double p = 0, q = 0;
  int T = 0;
  BitArray r = BitArray::zeros(graph->vertices());

  std::ifstream f_real;
  f_real.open(BENCHMARK_PATH + "/realization_num_" +
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
        p = stod(items[2]);
        break;
      case Q_LINE:
        q = stod(items[2]);
        break;
      case T_LINE:
        T = stoi(items[2]);
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

  Realization realization(p, q, T, r, BitArray::zeros(graph->vertices()));
  return std::unique_ptr<SourceDetectionParams>(
      new SourceDetectionParams(graph, realization, simulations));
}

std::string SourceDetectionParams::summary() const {
  return std::to_string(realization_.p()) + "_" +
         std::to_string(realization_.q()) + "_" +
         std::to_string(graph_->vertices());
}

}  // namespace cmplx
