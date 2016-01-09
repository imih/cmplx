#include "simulator.h"

#include "../common/idqueue.h"
#include "../common/ivector.cc"

using cmplx::common::IDqueue;
using cmplx::common::IGraph;
using cmplx::common::IVector;
using cmplx::common::BitArray;

namespace cmplx {
namespace simul {

void Simulator::NaiveSIROneStep(IGraph &graph, SirParams &sir_params) {
  BitArray &I = sir_params.infected();
  BitArray &S = sir_params.susceptible();
  BitArray &R = sir_params.recovered();
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
}

void Simulator::NaiveSIR(IGraph &graph, SirParams &sir_params) {
  while (sir_params.time_steps() < sir_params.maxT()) {
    NaiveSIROneStep(graph, sir_params);
  }
}

} // namespace simul
} // namespace cmplx
