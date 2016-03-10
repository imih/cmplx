#ifndef CMPLX_COMMON_RANDOM_H
#define CMPLX_COMMON_RANDOM_H

#include <syscall.h>
#include <unistd.h>

#include <random>

namespace cmplx {
namespace common {

class Random {
public:
  Random(long long seed) : generator_(seed), prob_distribution_(0, 1) {}
  /* Execute one experiment over Bernoulli random variable with probabilty p and
  * return the outcome.
  */
  bool eventDraw(double probabilty) { return prandReal01() <= probabilty; }

  double prandReal01() {
    // 64b mt19937
    return prob_distribution_(generator_);
  }

private:
  std::mt19937_64 generator_;
  std::uniform_real_distribution<double> prob_distribution_;
};
} // namespace cmplx
} // namespace common

#endif // CMPLX_COMMON_RANDOM_H
