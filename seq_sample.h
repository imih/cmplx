#ifndef CMPLX_SEQMC_STAT_H
#define CMPLX_SEQMC_STAT_H

#include "common/igraph.h"
#include "common/realization.h"
#include "common/bit_array.h"
#include "common/sir_params.h"
#include "simul/simulator.h"

#include <vector>

namespace cmplx {

class SeqSample {
 public:
  SeqSample(int v, const common::Realization& realization);
  SeqSample(const SeqSample& seqSample);

  ~SeqSample() = default;

  bool match(const common::Realization& realization) const;

  const common::BitArray& infected() const { return infected_; }
  const common::BitArray& recovered() const { return recovered_; }

  int t() const { return t_; }
  double w() const { return w_; }
  double pi() const { return pi_; }
  double g() const { return g_; }

  void update(const common::BitArray& new_infected,
              const common::BitArray& new_recovered, double newG, double newPi);

 private:
  common::BitArray realization() const { return infected_ | recovered_; }

  common::BitArray infected_;
  common::BitArray recovered_;
  int t_;
  double w_;   // w_t(x_t)
  double pi_;  // pi_t(x_t | x_{t - 1})
  double g_;   // g_t(x_t | x_{t - 1})
};
}

#endif  // CMPLX_SEQMC_STAT_H
