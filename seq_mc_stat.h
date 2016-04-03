#ifndef CMPLX_SEQMC_STAT_H
#define CMPLX_SEQMC_STAT_H

#include "common/realization.h"
#include "common/bit_array.h"
#include "common/sir_params.h"

#include <vector>

namespace cmplx {

class SeqSample {
 public:
  SeqSample(int v, const common::Realization& realization);
  SeqSample(const SeqSample& seqSample);

  common::BitArray infected() const { return sir_params_.infected(); }
  common::BitArray recovered() const { return sir_params_.recovered(); }

  common::SirParams& sir_params() { return sir_params_; }
  common::SirParams sir_params_const() const { return sir_params_; }

  int t() const { return t_; }
  std::vector<double> w() const { return w_; }

  void addW(double u) {
    double lastW = w_.back();
    w_.push_back(lastW * u);
    t_++;
  }

 private:
  common::SirParams sir_params_;
  int t_;
  std::vector<double> w_;
};

class SeqMCStat {
 public:
  SeqMCStat(int v, const common::Realization& realization)
      : realization_(realization), v_(v) {}

 private:
  SeqSample StartingSample();
  common::Realization realization_;
  int v_;
};
}

#endif  // CMPLX_SEQMC_STAT_H
