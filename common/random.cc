#include "random.h"

#include <unistd.h>
#include <sys/time.h>

namespace cmplx {
namespace common {

    Random::Random() : prob_distribution_(0, 1) {
      struct timeval tv;
      gettimeofday(&tv, 0);
      srand(tv.tv_usec * getpid());
      generator_.seed(tv.tv_usec * getpid());
    }

    bool Random::eventDraw(double probabilty) { return prandReal01() <= probabilty; }

    double Random::prandReal01() {
      // 64b mt19937
      return prob_distribution_(generator_);
    }
} // namespace cmplx
} // namespace common

