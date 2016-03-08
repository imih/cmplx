#ifndef CMPLX_COMMON_RANDOM_H
#define CMPLX_COMMON_RANDOM_H

#include <syscall.h>
#include <unistd.h>

#include <random>

namespace cmplx {
namespace common {

class Random {
public:
  Random(int64_t seed = time(NULL)) : seed_(seed) {}

  /* Execute one experiment over Bernoulli random variable with probabilty p and
  * return the outcome.
  */
  bool eventDraw(double probabilty) const {
    return prandReal01() <= probabilty;
  }

  double prandReal01() const {
    // 64b mt19937
    static thread_local std::mt19937_64 generator(seed_);
    static thread_local std::uniform_real_distribution<double>
        prob_distribution_(0, 1);
    return prob_distribution_(generator);
  }

private:
  int64_t seed_;
};
} // namespace cmplx
} // namespace common

#endif // CMPLX_COMMON_RANDOM_H
