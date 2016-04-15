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

  common::BitArray realization() const { return realization_; }

  bool match(const common::Realization& realization) const;

  int t() const { return t_; }
  double w() const { return w_; }
  double pi() const { return pi_; }
  double g() const { return g_; }

  void update(const common::BitArray new_realization, double newG,
              double newPi);

 private:
  common::BitArray realization_;
  int t_;
  double w_;   // w_t(x_t)
  double pi_;  // pi_t(x_t | x_{t - 1})
  double g_;   // g_t(x_t | x_{t - 1})
};

class SeqMCStat {
 public:
  SeqMCStat(int v, const common::Realization& target_realization,
            const common::IGraph& graph);

  void SISStep();

 private:
  common::BitArray draw_g_1(const common::BitArray& cur_realization);
  double g_1_cond(const common::BitArray& new_r, const common::BitArray& old_r);
  double Pi_1_cond(const common::BitArray& new_r,
                   const common::BitArray& old_r);

  void printvc2();

  std::vector<SeqSample> prev_samples_;
  std::vector<SeqSample> samples_;

  common::Realization target_realization_;
  std::vector<int> target_infected_idx_;
  int v_;
  cmplx::simul::Simulator simulator_;
};
}

#endif  // CMPLX_SEQMC_STAT_H
