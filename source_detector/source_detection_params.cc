#include "source_detection_params.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <random>

#include "../simul/simulator.h"
#include "../common/realization.h"

using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::Realization;
using cmplx::common::RealizationRead;

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
std::string BENCHMARK_PATH = "/home/imiholic";
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

  int simulations = 1e9;
  common::RealizationRead read(r, p, q, maxT);
  return std::unique_ptr<SourceDetectionParams>(
      new SourceDetectionParams(graph, read, simulations));
}

std::unique_ptr<SourceDetectionParams> SourceDetectionParams::LatticeCenter() {
  int lattice_size = 3;
  IGraph* graph = IGraph::UndirectedLattice({lattice_size, lattice_size});
  BitArray r = BitArray::ones(graph->vertices());
  double p = 0.2;
  double q = 0;
  int maxT = 2;
  common::RealizationRead read(r, p, q, maxT);
  return std::unique_ptr<SourceDetectionParams>(
      new SourceDetectionParams(graph, read, 1000000));
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
  int maxT = 5;
  int source_v = n * n / 2;
  BitArray zeros = BitArray::zeros(graph->vertices());

  cmplx::simul::Simulator simulator(graph);
  while (true) {
    BitArray infected = zeros;
    infected.set(source_v, true);
    BitArray susceptible = BitArray::ones(graph->vertices());
    Realization sir_params(p, q, maxT, susceptible, infected, zeros);
    simulator.NaiveSIR(sir_params);
    if (sir_params.realization().bitCount() > 1) {
      common::RealizationRead read(sir_params.realization(), p, q, maxT);
      return std::unique_ptr<SourceDetectionParams>(
          new SourceDetectionParams(graph, read, 1000000));
    }
  }
}

std::unique_ptr<SourceDetectionParams> SourceDetectionParams::ParamsFromGridISS(
    double p, double q, int n) {
  IGraph* graph = IGraph::UndirectedLattice({n, n});
  int maxT = 5;
  int source_v = n * n / 2;
  BitArray zeros = BitArray::zeros(graph->vertices());

  cmplx::simul::Simulator simulator(graph);
  while (true) {
    BitArray infected = zeros;
    infected.set(source_v, true);
    BitArray s = BitArray::ones(graph->vertices());
    Realization sir_params(p, q, maxT, s, infected, zeros);
    simulator.NaiveISS(sir_params);
    if (sir_params.realization().bitCount() > 1) {
      common::RealizationRead read(sir_params.realization(), p, q, maxT);
      return std::unique_ptr<SourceDetectionParams>(
          new SourceDetectionParams(graph, read, 1000000));
    }
  }
}

std::unique_ptr<SourceDetectionParams> SourceDetectionParams::ParamsFromGML(
    const std::string& file_name, int source_node, double p, double q) {
  IGraph* graph = IGraph::GraphFromGML(file_name);
  int maxT = 5;
  BitArray zeros = BitArray::zeros(graph->vertices());

  cmplx::simul::Simulator simulator(graph);
  while (true) {
    BitArray infected = zeros;
    infected.set(source_node, true);
    BitArray s = BitArray::ones(graph->vertices());
    Realization sir_params(p, q, maxT, s, infected, zeros);
    simulator.NaiveSIR(sir_params);
    if (sir_params.realization().bitCount() > 1) {
      common::RealizationRead read(sir_params.realization(), p, q, maxT);
      SourceDetectionParams* params =
          new SourceDetectionParams(graph, read, 1000000);
      params->setSourceId(source_node);
      return std::unique_ptr<SourceDetectionParams>(params);
    }
  }
}

std::unique_ptr<SourceDetectionParams> SourceDetectionParams::BenchmarkParams(
    int realization_no) {
  while (realization_no > 160) realization_no -= 160;

  IGraph* graph =
      IGraph::GraphFromGML(BENCHMARK_PATH + "/network/lattice30.gml");
  double p = 0, q = 0;
  int T = 0;
  BitArray r = BitArray::zeros(graph->vertices());

  std::ifstream f_real;
  std::string filename = BENCHMARK_PATH + "/realizations/realization_" +
                         std::to_string(realization_no) + ".txt";
  f_real.open(filename);
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
  common::RealizationRead read(r, p, q, T);
  return std::unique_ptr<SourceDetectionParams>(
      new SourceDetectionParams(graph, read, 100000));
}

std::string SourceDetectionParams::summary() const {
  std::string s = std::to_string(realization_.p()) + "_" +
                  std::to_string(realization_.q()) + "_" +
                  std::to_string(graph_->vertices());
  return s;
}

}  // namespace cmplx
