#ifndef SOURCE_DETECTION_PARAL_H
#define SOURCE_DETECTION_PARAL_H

#include "source_detection_params.h"
#include "source_detector.h"
#include <vector>

namespace cmplx {

void DirectMCSimulParalConv(const cmplx::SourceDetectionParams& params,
                            ModelType modelType = ModelType::SIR);
void DirectMCSimulParal(const cmplx::SourceDetectionParams& params,
                        ModelType modelType = ModelType::SIR);

std::vector<double> SoftMarginParalConvMaster(
    const cmplx::SourceDetectionParams& params, bool end = true);
void SoftMarginParalConv(const cmplx::SourceDetectionParams& params,
                         ModelType modeltype = ModelType::SIR);
void SoftMarginParal(const cmplx::SourceDetectionParams& params,
                     ModelType modelType = ModelType::SIR);

void GenerateSoftMarginDistributions(const cmplx::SourceDetectionParams& params,
                                     int distributions,
                                     ModelType modelType = ModelType::SIR);

void GenerateSeqMonteCarloDistributions(
    const cmplx::SourceDetectionParams& params, int distributions);
std::vector<double> SeqMonteCarloSimulParalMaster(
    const SourceDetectionParams &params, bool end, bool print);
void SeqMonteCarloSimulParalWorker(const SourceDetectionParams &params);


}  // namespace cmplx
#endif  // SOURCE_DETECTION_PARAL_H
