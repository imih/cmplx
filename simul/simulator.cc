#include "simulator.h"

#include "../common/idqueue.h"
#include "../common/ivector.cc"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

using cmplx::common::IDqueue;
using cmplx::common::IGraph;
using cmplx::common::IVector;
using cmplx::common::BitArray;
using cmplx::common::SirParams;

typedef std::numeric_limits<double> dbl;

namespace cmplx {
namespace simul {

void Simulator::NaiveSIROneStep(const IGraph &graph, SirParams &sir_params) {
  BitArray I = sir_params.infected();
  BitArray S = sir_params.susceptible();
  BitArray R = sir_params.recovered();
  IDqueue &infected_q = sir_params.infected_q();
  int batch_size = I.bits_num();

  while (!infected_q.empty() && batch_size) {
    int u = infected_q.pop();
    batch_size--;
    const IVector<int> &adj_list_u = graph.adj_list(u);
    int adj_list_size = adj_list_u.size();
    for (int idx = 0; idx < adj_list_size; ++idx) {
      int v = adj_list_u[idx];
      if (S.bit(v) && sir_params.drawP()) {
        S.set(v, 0);
        I.set(v, 1);
        infected_q.push(v);
      }
    }
    if (sir_params.drawQ()) {
      I.set(u, 0);
      R.set(u, 1);
    }
  }
  sir_params.incrTime();
  sir_params.set_infected(I);
  sir_params.set_susceptible(S);
  sir_params.set_recovered(R);
}

void Simulator::NaiveSIR(const IGraph &graph, SirParams &sir_params) {
  while (sir_params.time_steps() < sir_params.maxT()) {
    NaiveSIROneStep(graph, sir_params);
  }
}

namespace {
double *bin_coef_cache;
void build_bin_coef_cache(int n) {
  bin_coef_cache = new double[(n + 1) * (n + 1)];
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n; ++j) {
      double &ref = bin_coef_cache[(n + 1) * i + j];
      if (!j || i == j)
        ref = 1;
      else if (i < j)
        ref = 0;
      else if (i && j) {
        ref = bin_coef_cache[(n + 1) * (i - 1) + j - 1] +
              bin_coef_cache[(n + 1) * (i - 1) + j];
      } else
        ref = 0;
    }
}

void clear_bin_coef_cache() { delete bin_coef_cache; }
double Ck(int n, int k) { return bin_coef_cache[(n + 1) * n + k]; }

double P(int n, int k, double p, double q) {
  double ret = 0;
  for (int l = 0; l <= k; ++k) {
    double cur = l % 2 ? -1 : 1;
    double powP = pow(1 - p, n - k + l);
    cur *= Ck(k, l) * powP;
    cur /= (1 - (1 - q) * powP);
    ret += cur;
  }
  return ret * q * Ck(n, k);
}

} // anonymous namespace for CDF calculations

void Simulator::calcCummulativeInfecting(int n, const std::string &file_name) {
  build_bin_coef_cache(n);

  std::ofstream ost;
  ost.open(file_name.c_str());
  ost.precision(dbl::max_digits10);
  ost << n << std::endl;
  for (double p = 0; p <= 1; p += 0.1) {
    for (double q = 0; q <= 1; q += 0.1) {
      ost << p << " " << q << std::endl;
      double prev = 0;
      for (int k = 0; k <= n; ++k) {
        double pK = P(n, k, p, q);
        prev += pK;
        ost << prev << " ";
      }
      ost << std::endl;
    }
  }
  ost.close();
  clear_bin_coef_cache();
}

} // namespace simul
} // namespace cmplx
