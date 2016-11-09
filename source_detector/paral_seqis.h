#ifndef PARAL_SEQIS_H
#define PARAL_SEQIS_H

#include "source_detection_params.h"
#include "common_paral.h"

#include <vector>
#include <string>

using std::vector;

namespace cmplx {
class ParalSeqIS : public CommonTraits {
 public:
  ParalSeqIS()
      : CommonTraits({(int)1e2, (int)1e3, (int)1e4, (int)1e5, (int)1e6},
                     "SEQSoftbench_SBS_",
                     {(int)1e4,      2 * (int)1e4,  4 * (int)1e4,
                      8 * (int)1e4,  10 * (int)1e4, 20 * (int)1e4,
                      40 * (int)1e4, 80 * (int)1e4, 100 * (int)1e4}) {}

  ~ParalSeqIS() = default;
};
}

#endif  // PARAL_SEQIS_H
