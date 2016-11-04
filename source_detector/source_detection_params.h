#ifndef SOURCE_DETECTION_PARAMS_H
#define SOURCE_DETECTION_PARAMS_H

#include "../common/igraph.h"
#include "../common/realization.h"
#include "../common/bit_array.h"
#include <string>
#include <memory>

namespace cmplx {
class SourceDetectionParams {
 public:
  // TODO FIX!
  static std::unique_ptr<SourceDetectionParams> SupFig2Params();
  static std::unique_ptr<SourceDetectionParams> LatticeCenter();
  static std::unique_ptr<SourceDetectionParams> BenchmarkParams(
      int realization_no);
  static std::unique_ptr<SourceDetectionParams> ParamsFromGrid(double p,
                                                               double q, int n);
  static std::unique_ptr<SourceDetectionParams> ParamsFromGridISS(double p,
                                                                  double q,
                                                                  int n);

  static std::unique_ptr<SourceDetectionParams> ParamsFromGML(
      const std::string &file_name, int source_node, double p, double q);

  ~SourceDetectionParams() = default;

  const std::unique_ptr<common::IGraph> &graph() const { return graph_; }

  const common::RealizationRead &realization() const { return realization_; }

  long long simulations() const { return simulations_; }

  double a() const { return a_; }

  int sourceID() { return source_id_; }

  void setRealization(const common::BitArray &r) {
    realization_.setRealization(r);
  }

  void setSimulations(long long simulations) { simulations_ = simulations; }

  void setA(double a) { a_ = a; }

  std::string summary() const;
  void setSourceId(int source_id) { source_id_ = source_id; }

 private:
  SourceDetectionParams(common::IGraph *graph,
                        const common::RealizationRead &realization,
                        long long simulations, double a = pow(2, -5))
      : graph_(graph),
        realization_(realization),
        a_(a),
        simulations_(simulations) {}

  std::unique_ptr<common::IGraph> graph_;
  common::RealizationRead realization_;
  double a_;
  long long simulations_;
  int source_id_;
};

}  // namespace cmplx

#endif  // SOURCE_DETECTION_PARAMS_H
