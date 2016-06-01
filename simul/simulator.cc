#include "simulator.h"

#include "../common/idqueue.h"
#include "../common/ivector.h"
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
  assert(S.bitCount() == S.bits_num() - 1);
  IDqueue q(I);
  assert(q.size() == 1);

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
      if(!I.bit(current_node)) {
        puts("DOGODILO SE");
        continue;
      }

      const IVector<int> &neis = graph_->adj_list(current_node);
      for (int i = 0; i < neis.size(); ++i) {
        int current_neigh = neis[i];
        if (S.bit(current_neigh)) {
          if (eventDraw(sir_params.p())) {
            q.push(current_neigh);
            I.set(current_neigh, true);
            S.set(current_neigh, false);
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

void Simulator::NaiveSIROneStep(common::SirParams &sir_params) {
  puts("DEPRECATED");
  exit(1); 
  BitArray I = sir_params.infected();
  BitArray R = sir_params.recovered();
  BitArray S = sir_params.susceptible();
  double p = sir_params.p();
  double q = sir_params.q();

  IDqueue queue(I);

  int delta_nodes_pop = queue.size();
  int p0 = 0, p1 = 0, q0 = 0, q1 = 0;

  while (queue.size()) {
    if (delta_nodes_pop == 0) {
      break;
    }
    long int current_node = queue.pop();
    delta_nodes_pop--;

    const IVector<int> &neis = graph_->adj_list(current_node);
    for (int i = 0; i < neis.size(); ++i) {
      int current_neigh = neis[i];

      if (S.bit(current_neigh)) {
        if (eventDraw(p)) {
          p1++;
          queue.push(current_neigh);
          I.set(current_neigh, true);
          S.set(current_neigh, false);
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
  sir_params.set_susceptible(S);
}

bool Simulator::NaiveISS(common::SirParams &sir_params, bool prunning,
                         const common::BitArray &allowed_nodes) {
  // p . lambda
  // q . alfa
  // S . ignorant
  // I . spreader
  // R . stifler
  BitArray I = sir_params.infected();
  BitArray S = sir_params.susceptible();
  BitArray R = sir_params.recovered();

  assert(I.bitCount() == 1);
  assert(R.bitCount() == 0);
  IDqueue q(I);  // sir_params.population_size());
  assert(q.size() == 1);
  assert(S.bitCount() == S.bits_num() - 1);

  int dis_time = 1;
  int delta_nodes_pop = q.size();
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
      if (!I.bit(current_node)) {
        continue;
      }

      const IVector<int> &neis = graph_->adj_list(current_node);
      bool became_stifler = false;
      for (int i = 0; i < neis.size(); ++i) {
        int current_neigh = neis[i];
        if (S.bit(current_neigh)) {
          // spreader meets ignorant
          if (eventDraw(sir_params.p())) {
            S.set(current_neigh, false);
            I.set(current_neigh, true);
            q.push(current_neigh);
            if (prunning && allowed_nodes.bit(current_neigh) == 0) {
              prunned = true;
              break;
            }
          }
        } else {
          if (eventDraw(sir_params.q())) {
            became_stifler = true;
            I.set(current_node, false);
            R.set(current_node, true);
            //break;
          }
        }
      }
      if (!became_stifler) {
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
