#ifndef SOURCE_DETECTION_PARAL_H
#define SOURCE_DETECTION_PARAL_H

#include "source_detection_params.h"
#include <vector>

namespace cmplx {

void DirectMCSimulParalConv(const cmplx::SourceDetectionParams& params);
void DirectMCSimulParal(const cmplx::SourceDetectionParams& params);

std::vector<double> SoftMarginParalConvMaster(
    const cmplx::SourceDetectionParams& params, bool end = true);
void SoftMarginParalConv(const cmplx::SourceDetectionParams& params);
void SoftMarginParal(const cmplx::SourceDetectionParams& params);

void GenerateSoftMarginDistributions(const cmplx::SourceDetectionParams& params,
                                    int distributions);

// estimates the posterior probabilty of a full match
void FEstimatorParal(const cmplx::SourceDetectionParams& params);

}  // namespace cmplx
#endif  // SOURCE_DETECTION_PARAL_H
