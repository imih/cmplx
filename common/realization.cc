#include "realization.h"

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
}  // namespace

namespace cmplx {
namespace common {
Realization::Realization(double p, double q, int maxT,
                         const BitArray& susceptible, const BitArray& infected,
                         const BitArray& recovered)
    : p_(p),
      q_(q),
      T_(maxT),
      susceptible_(susceptible),
      infected_(infected),
      recovered_(recovered),
      realization_(infected | recovered) {}

void Realization::set_infected(const BitArray& infected) {
  infected_ = infected;
  calcRealization();
}

void Realization::set_recovered(const BitArray& recovered) {
  recovered_ = recovered;
  calcRealization();
}

void Realization::update(const BitArray& infected, const BitArray& recovered) {
  infected_ = infected;
  recovered_ = recovered;
  calcRealization();
}

void Realization::printForLattice(int n) const {
  string I = infected_.to_string();
  string S = susceptible_.to_string();
  string R = recovered_.to_string();
  std::cout << "I\n" << split(I, n) << " S\n" << split(S, n) << " R\n"
            << split(R, n) << std::endl;
}

}  // namespace common
}  // namespace cmplx
