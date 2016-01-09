#include "sir_params.h"

namespace cmplx {
namespace common {
SirParams::SirParams(double p, double q, int T, const BitArray &infected,
                     const BitArray &susceptible)
    : random_(), infected_q_(infected), time_steps_(0), p_(p), q_(q), T_(T),
      infected_(infected), susceptible_(susceptible),
      recovered_(infected.bits_num()) {}

} // namespace common
} // namespace cmplx

