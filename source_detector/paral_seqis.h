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
      : CommonTraits({100, 200, 1000, 2000, 
                     (int)1e4, 2 * (int)1e4,  
                      (int)1e5, 2 * (int)1e5,
                      (int)1e6, 2 * (int) 1e6},
                     "SEQbench_",
                     {(int)1e4, 2 * (int)1e4,  4 * (int)1e4,
                      (int)1e5, 2 * (int)1e5,
                      4 * (int)1e5 , (int)1e6, 2 * (int) 1e6, 4 * (int) 1e6}) {}

  ~ParalSeqIS() = default;
};

class ParalSoftSeqIS : public CommonTraits {
 public:
  ParalSoftSeqIS()
      : CommonTraits({100, 200, 1000, 2000, 
                     (int)1e4, 2 * (int)1e4,  
                      (int)1e5, 2 * (int)1e5,
                      (int)1e6, 2 * (int) 1e6},
                     "SEQSoftbench_",
                     {(int)1e4, 2 * (int)1e4,  4 * (int)1e4,
                      (int)1e5, 2 * (int)1e5,
                      4 * (int)1e5 , (int)1e6, 2 * (int) 1e6, 4 * (int) 1e6}) {}

  ~ParalSoftSeqIS() = default;
};
}

#endif  // PARAL_SEQIS_H
