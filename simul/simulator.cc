#include "simulator.h"

#include "../common/idqueue.h"

using cmplx::common::IDqueue;
using cmplx::common::IGraph;
using cmplx::common::IVector;
using cmplx::common::BitArray;

namespace cmplx {
namespace simul {
SimulationStats Simulator::NaiveSIR(IGraph &graph, SirParams &sir_params) {
  SimulationStats ss(sir_params);
  BitArray &I = ss.infected();
  BitArray &S = ss.susceptible();
  BitArray &R = ss.recovered();
  IDqueue infected_q(I);

  while (!infected_q.empty()) {
    //    dequeue(u, I);
    int u = infected_q.pop();

    const IVector<int> &adj_list_u = graph.adj_list(u);
    int adj_list_size = adj_list_u.size();
    for (int idx = 0; idx < adj_list_size; ++idx) {
      int v = adj_list_u[idx];
      if (S.bit(v)) {
        /*           let the transmission of infection u->v occur with
        * probabilty p
        *           if u->v occurs:
        *              update S(v) and I(v)
        *              enqueue(v,I)
        *              */
      }
    }
    /*
     *      update state u from infected to recovered with probabilty q
     *      if u is not recovered
     *         enqueue(u, I)
     */
  }
  return ss;
}

} // namespace simul
} // namespace cmplx
