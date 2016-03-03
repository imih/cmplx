#ifndef CMPLX_COMMON_RANDOM_H
#define CMPLX_COMMON_RANDOM_H

#include <syscall.h>
#include <unistd.h>

#include <random>

namespace cmplx {
namespace common {
class Random {
public:
  Random() {}

  /* Execute one experiment over Bernoulli random variable with probabilty p and
  * return the outcome.
  */
  static bool eventDraw(double probabilty) {
    static thread_local std::mt19937 generator(randomInt32());
    static std::uniform_real_distribution<double> prob_distribution_(0, 1);
    // 32b mt19937
    return prob_distribution_(generator) <= probabilty;
  }

  static int64_t randomInt32() {
    int64_t r;
    // from /dev/urandom
    syscall(SYS_getrandom, (void *)&r, sizeof r, 0);
    return r;
  }
};
} // namespace cmplx
} // namespace common

#endif // CMPLX_COMMON_RANDOM_H
