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

bool Simulator::NaiveSIR(SirParams &sir_params, bool prunning,
                         const BitArray &allowed_nodes) {
  BitArray I = sir_params.infected();
  BitArray S = sir_params.susceptible();
  BitArray R = sir_params.recovered();

  assert(I.bitCount() == 1);
  IDqueue q(sir_params.population_size());
  q.insertMarked(I);

  int dis_time = 1;
  int delta_nodes_pop = q.size();
  long int num_inf_nodes = 1;  // source
  bool prunned = false;

  while (q.size() && !prunned) {
    if (delta_nodes_pop == 0) {
      dis_time++;
      delta_nodes_pop = q.size();
    }

    if (dis_time <= sir_params.maxT()) {
      long int current_node;
      current_node = q.pop();
      delta_nodes_pop--;

      const IVector<int> &neis = graph_.adj_list(current_node);
      for (int i = 0; i < neis.size(); ++i) {
        int current_neigh = neis[i];

        if (I.bit(current_neigh) == 0 && R.bit(current_node) == 0) {
          if (eventDraw(sir_params.p())) {
            q.push(current_neigh);
            I.set(current_neigh, true);
            num_inf_nodes++;
            // indentity prunning
            if (prunning && allowed_nodes.bit(current_neigh) == 0) {
              prunned = true;
              break;
            }
          }
        }
      }

      if (eventDraw(sir_params.q())) {
        R.set(current_node, true);
        I.set(current_node, false);
      } else {
        q.push(current_node);
      }
    } else {
      break;
    }
  }

  sir_params.set_infected(I);
  sir_params.set_susceptible(S);
  sir_params.set_recovered(R);
  return prunned;
}

double Simulator::NaiveSIROneStep(common::SirParams &sir_params) {
  BitArray I = sir_params.infected();
  BitArray R = sir_params.recovered();
  double p = sir_params.p();
  double q = sir_params.q();

  IDqueue queue(sir_params.population_size());
  queue.insertMarked(I);
  int delta_nodes_pop = queue.size();
  int p0 = 0, p1 = 0, q0 = 0, q1 = 0;

  while (queue.size()) {
    if (delta_nodes_pop == 0) {
      break;
    }
    long int current_node = queue.pop();
    delta_nodes_pop--;

    const IVector<int> &neis = graph_.adj_list(current_node);
    for (int i = 0; i < neis.size(); ++i) {
      int current_neigh = neis[i];

      if (I.bit(current_neigh) == 0 && R.bit(current_node) == 0) {
        if (eventDraw(p)) {
          p1++;
          queue.push(current_neigh);
          I.set(current_neigh, true);
          // indentity prunning
        } else
          p0++;
      }
    }

    if (eventDraw(q)) {
      q1++;
      R.set(current_node, true);
      I.set(current_node, false);
    } else {
      q0++;
      queue.push(current_node);
    }
  }

  sir_params.set_infected(I);
  sir_params.set_recovered(R);
  return pow(1 - p, p0) * pow(p, p1) * pow(1 - q, q0) * pow(q, q1);
}

/*
BitArray I = sir_params.infected();
BitArray S = sir_params.susceptible();
BitArray R = sir_params.recovered();
IDqueue infected_q(sir_params.population_size());
infected_q.clear();
infected_q.insertMarked(I);

int t = 1;
int batch_size = infected_q.size();
bool prunned = false;

while (!infected_q.empty() && !prunned) {
  if (batch_size == 0) {
    t++;
    batch_size = infected_q.size();
  }
  if (t > sir_params.maxT()) {
    break;
  }

  int u = infected_q.pop();
  batch_size--;
  const IVector<int> &adj_list_u = graph_.adj_list(u);

  int adj_list_size = adj_list_u.size();
  for (int idx = 0; idx < adj_list_size; ++idx) {
    int v = adj_list_u[idx];
    if (S.bit(v) && draw(sir_params.p())) {
      S.set(v, 0);
      I.set(v, 1);
      infected_q.push(v);
      if (prunning && !allowed_nodes.bit(v)) {
        prunned = true;
        break;
      }
    }
  }
  if (draw(sir_params.q())) {
    I.set(u, 0);
    R.set(u, 1);
  } else {
    infected_q.push(u);
  }
}

sir_params.set_infected(I);
sir_params.set_susceptible(S);
sir_params.set_recovered(R);
infected_q.clear();
return prunned;
*/

/*
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
*/

/*
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
*/

}  // namespace simul
}  // namespace cmplx
