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
      : CommonTraits({100, 200, 1000, 2000, 
                     (int)1e4, 2 * (int)1e4,  
                      (int)1e5, 2 * (int)1e5,
                      (int)1e6, 2 * (int) 1e6},
            "SMbench_",
            {(int)1e4,       2 * (int)1e4,   4 * (int)1e4,   (int)1e5,
             2 * (int)1e5,  4 * (int)1e5,  (int)1e6, 2 * (int)1e6,
             4 * (int)1e6, (int)1e7}) {}

  ~ParalSoftMC() = default;
};
}

#endif  // PARAL_SOFTMC_H
