#include "simulator.h"

#include "../common/idqueue.h"
#include "../common/ivector.h"
#include "../common/ivector.cc"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <map>

using cmplx::common::IDqueue;
using cmplx::common::IGraph;
using cmplx::common::IVector;
using cmplx::common::BitArray;
using cmplx::common::SirParams;

typedef std::numeric_limits<double> dbl;

namespace cmplx {
namespace simul {

Simulator::Simulator(const common::IGraph *graph)
    : graph_(graph), prob_distribution_(0, 1) {
  struct timeval tv;
  gettimeofday(&tv, 0);
  srand(tv.tv_usec * getpid());
  generator_.seed(tv.tv_usec * getpid());
}

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
  IDqueue q(I);
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

}  // namespace simul
}  // namespace cmplx
