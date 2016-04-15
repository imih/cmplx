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

class SeqMCStat {
 public:
  SeqMCStat(int v, const common::Realization& target_realization,
            const common::IGraph& graph);

  double run() {
    for (int t = 0; t < target_realization_.maxT(); ++t) {
      SISStep();
    }

    double pos_P = 0;
    for (const SeqSample& sample : samples_) {
      if (sample.match(target_realization_)) {
        pos_P += sample.w();
      }
    }
    printf("%.10lf\n", pos_P);
    return pos_P;
  }

 private:
  void SISStep();

  std::pair<common::BitArray, common::BitArray> draw_g_1(
      const common::BitArray& cur_inf, const common::BitArray& cur_rec);
  double g_1_cond(const common::BitArray& new_i,
                  const common::BitArray& new_rec,
                  const common::BitArray& old_i,
                  const common::BitArray& old_rec);
  double Pi_1_cond(const common::BitArray& new_i,
                   const common::BitArray& new_rec,
                   const common::BitArray& old_i,
                   const common::BitArray& old_rec);

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
