#ifndef PARAL_DIRECTMC_H
#define PARAL_DIRECTMC_H

#include "source_detection_params.h"
#include "common_paral.h"

#include <vector>
#include <string>

using std::vector;

namespace cmplx {
class ParalDirectMC : public CommonTraits {
 public:
  ParalDirectMC()
      : CommonTraits(
            {(int)1e4, (int)1e5, (int)1e6, (int)1e7, (int)1e8, (int)1e9},
            "DMCbench_",
            {(int)1e6,     2 * (int)1e6, 4 * (int)1e6, (int)1e7,
             2 * (int)1e7, 4 * (int)1e7, (int)1e8,     2 * (int)1e8,
             4 * (int)1e8, (int)1e9}) {}
  ~ParalDirectMC() = default;
};
}

#endif  // PARAL_DIRECTMC_H
