#include "./source_detector.h"
#include "common/bit_array.h"
#include "common/igraph.h"
#include "common/realization.h"
#include "common/sir_params.h"
#include "source_detection_params.h"

#include <iostream>
#include <vector>
#include <string>

#include <igraph/igraph.h>

using cmplx::SourceDetector;
using cmplx::common::IGraph;
using cmplx::common::BitArray;
using cmplx::common::SirParams;
using cmplx::common::Realization;
using cmplx::SourceDetectionParams;

/*
void generateBarabasi(int no, int poc = 0) {
  for (int t = poc; t < poc + no; ++t) {
    printf("%d\n", t);
    std::string file_name = "barabasi2_900_" + std::to_string(t) + ".gml";
    IGraph* g = IGraph::BarabasiAlbert(900);
    g->writeGML(file_name);
    std::string info_name = "barabasi2_900_" + std::to_string(t) + ".info";
    FILE* f = fopen(info_name.c_str(), "w");
    fprintf(f, "id,deg,kcore,clos,betw,eigcentr\n");
    for (int i = 0; i < g->vertices(); ++i) {
      fprintf(f, "%d,%d,%d,%.10lf,%.10lf,%.10lf\n", i, g->deg(i), g->kCore(i),
              g->closeness(i), g->betweenness(i), g->eigenvector_centrality(i));
    }
    fclose(f);
    delete (g);
  }
}

void generateErdosRenyi(int no, int poc = 0) {
  for (int t = poc; t < poc + no; ++t) {
    int gs = 0;
    while (true) {
      printf("gs: %d\n", gs);
      printf("%d\n", t);
      IGraph* g = IGraph::ErdosRenyi(900, 0.01);
      if (g->is_connected()) {
        std::string file_name =
            "erdos_renyi_900_0.01_" + std::to_string(t) + ".gml";
        g->writeGML(file_name);
        std::string info_name =
            "erdos_renyi_900_0.01_" + std::to_string(t) + ".info";
        FILE* f = fopen(info_name.c_str(), "w");
        fprintf(f, "id,deg,kcore,clos,betw,eigcentr\n");
        for (int i = 0; i < g->vertices(); ++i) {
          fprintf(f, "%d,%d,%d,%.10lf,%.10lf,%.10lf\n", i, g->deg(i),
                  g->kCore(i), g->closeness(i), g->betweenness(i),
                  g->eigenvector_centrality(i));
        }
        fclose(f);
        delete g;
        break;
      }
      delete g;
      gs++;
    }
  }
}
*/

int main() {
  // generateBarabasi(50, 0);
  // generateErdosRenyi(50, 0);

  // fclose(f);

  auto params = cmplx::SourceDetectionParams::BenchmarkParams(1);
  /*
  for (double p : probs) {
    printf("%.10lf\n", p);
  }
  */
  return 0;
}
