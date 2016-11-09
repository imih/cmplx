#ifndef PARAL_SOFTMC_H
#define PARAL_SOFTMC_H

#include "source_detection_params.h"
#include "common_master.h"

#include <vector>
#include <string>

using std::vector;

namespace cmplx {
class ParalSoftMC : public CommonTraits {
 public:
  ParalSoftMC()
      : CommonTraits(
            {(int)1e2, (int)1e3, (int)1e4, (int)1e5, (int)1e6}, "SMbench_SBS_",
            {(int)1e4,       2 * (int)1e4,   4 * (int)1e4,   10 * (int)1e4,
             20 * (int)1e4,  40 * (int)1e4,  100 * (int)1e4, 200 * (int)1e4,
             400 * (int)1e4, 1000 * (int)1e4}) {}

  ~ParalSoftMC() = default;
};
}

#endif  // PARAL_SOFTMC_H
