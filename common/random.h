#ifndef CMPLX_COMMON_RANDOM_H
#define CMPLX_COMMON_RANDOM_H

#include <random>

namespace cmplx {
namespace common {

class Random {
public:
  Random();

  /* Execute one experiment over Bernoulli random variable with probabilty p and
  * return the outcome.
  */
  bool eventDraw(double probabilty);

  double prandReal01();

private:
  std::mt19937_64 generator_;
  std::uniform_real_distribution<double> prob_distribution_;
};
} // namespace cmplx
} // namespace common

#endif // CMPLX_COMMON_RANDOM_H
