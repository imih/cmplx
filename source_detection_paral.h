#ifndef SOURCE_DETECTION_PARAL_H
#define SOURCE_DETECTION_PARAL_H

#include "source_detection_params.h"

namespace cmplx {
void DirectMCSimulParal(cmplx::SourceDetectionParams params);

void SoftMarginParal(cmplx::SourceDetectionParams params);

// estimates the posterior probabilty of a full match
void FEstimatorParal(cmplx::SourceDetectionParams params);

}  // namespace cmplx
#endif  // SOURCE_DETECTION_PARAL_H
