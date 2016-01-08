#include "simulator.h"

using cmplx::common::IGraph;
using cmplx::common::IVector;

namespace cmplx {
SimulationStats Simulator::NaiveSIRStep(const IGraph& graph,
                                   const SirParams& sir_params) {
  /*
   * I -- queue of infected nodes
   *  while(!I.empty()) {
   *     dequeue(u, I);
   *     for each contact v of node u:
   *        if(S(v)):
   *           let the transmission of infection u->v occur with probabilty p
   *           if u->v occurs:
   *              update S(v) and R(v)
   *              enqueue(v,I)
   *      update state u from infected to recovered with probabilty q
   *      if u is not recovered 
   *         enqueue(u, I)
   */

  // TODO
  SimulationStats ss;
  return ss;
}

/*
SimulationStats simulateSIRNaive(const IGraph& graph,
                                        const IVector<int>& infected_nodes,
                                        const SirParams& sir_params) {
  //TODO
  SimulationStats ss;
  return ss;
}
*/
}  // namespace cmplx
