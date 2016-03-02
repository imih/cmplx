#ifndef CMPLX_COMMON_RANDOM_H
#define CMPLX_COMMON_RANDOM_H

#include <random>

namespace cmplx {
namespace common {
class Random {
public:
  Random() : prob_distribution_(0, 1) {}

  /* Execute one experiment over Bernoulli random variable with probabilty p and
  * return the outcome.
  */
  bool eventDraw(double probabilty) {
    static thread_local std::mt19937 generator;
    return prob_distribution_(generator) <= probabilty;
  }

private:
  std::uniform_real_distribution<double> prob_distribution_;
};
} // namespace cmplx
} // namespace common

#endif // CMPLX_COMMON_RANDOM_H
