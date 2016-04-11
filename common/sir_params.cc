#include "sir_params.h"

#include <string>
#include <sstream>

using std::string;

namespace {
  string split(string s, int n) {
    std::ostringstream oss;
    for (int i = 0; i < (int)s.size(); i += n) {
      oss << s.substr(i, n) << std::endl;
    }
    return oss.str();
  }
} // namespace

namespace cmplx {
  namespace common {
    SirParams::SirParams(double p, double q, int T, const BitArray &infected,
        const BitArray &susceptible)
      : p_(p), q_(q),
      T_(T), infected_(new BitArray(infected)), susceptible_(new BitArray(susceptible)),
      recovered_(new BitArray(infected.bits_num())) {}

    void SirParams::print() const {
      std::cout << "I " << infected_ << " S " << susceptible_ << " R " << recovered_
        << std::endl;
    }

    void SirParams::printForLattice(int n) const {
      string I = infected_->to_string();
      string S = susceptible_->to_string();
      string R = recovered_->to_string();
      std::cout << "I\n" << split(I, n) << " S\n" << split(S, n) << " R\n"
        << split(R, n) << std::endl;
    }

  } // namespace common
} // namespace cmplx
