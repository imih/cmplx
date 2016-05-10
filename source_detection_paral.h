#ifndef SOURCE_DETECTION_PARAL_H
#define SOURCE_DETECTION_PARAL_H

#include "source_detection_params.h"
#include "source_detector.h"
#include <vector>

namespace cmplx {

void DirectMCSimulParalConv(SourceDetectionParams* params,
                            ModelType modelType = ModelType::SIR);
void DirectMCSimulParal(const SourceDetectionParams* params,
                        ModelType modelType = ModelType::SIR);

void DirectMCBenchmark(SourceDetectionParams* params, int benchmark_no);


std::vector<double> SoftMarginParalConvMaster(SourceDetectionParams* params,
                                              bool end = true);
void SoftMarginParalConv(SourceDetectionParams* params,
                         ModelType modeltype = ModelType::SIR);
void SoftMarginParal(const SourceDetectionParams* params,
                     ModelType modelType = ModelType::SIR);

void GenerateSoftMarginDistributions(SourceDetectionParams* params,
                                     int distributions,
                                     ModelType modelType = ModelType::SIR);

void SoftMarginBenchmarkConv(SourceDetectionParams* params, int benchmark_no,
                             ModelType model_type = ModelType::SIR);

void GenerateSeqMonteCarloDistributions(SourceDetectionParams* params,
                                        int distributions);
void SeqMonteCarloBenchmark(SourceDetectionParams *params,
                                        int benchmark_no);
std::vector<double> SeqMonteCarloSimulParalMaster(
    const SourceDetectionParams* params, bool end, bool print);
void SeqMonteCarloSimulParalWorker(const SourceDetectionParams& params);

}  // namespace cmplx
#endif  // SOURCE_DETECTION_PARAL_H
