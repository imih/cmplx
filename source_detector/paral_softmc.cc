#include "paral_softmc.h"

#include <vector>
#include <string>
#include <algorithm>

namespace cmplx {

namespace {
double dabs(double x) {
  if (x < 0) return x * -1;
  return x;
}
}  // namespace

/*
vector<double> ParalSoftMC::softConvMaster(
    cmplx::SourceDetectionParams* params) {
  std::vector<int> sims = {(int)1e4,       2 * (int)1e4,   4 * (int)1e4,
                           10 * (int)1e4,  20 * (int)1e4,  40 * (int)1e4,
                           100 * (int)1e4, 200 * (int)1e4, 400 * (int)1e4,
                           800 * (int)1e4};

  double c = 0.05;
  const int MAXA = 9;
  int s0 = sims[0];
  printf("s0: %d\n", s0);
  vector<double> a(MAXA + 1, 0);
  for (int i = 3; i <= MAXA; ++i) {
    a[i] = 1.0 / (double)(1 << i);
  }
  vector<double> p0[MAXA + 1];
  vector<double> pMAP0(MAXA + 1, 0);
  for (int i = 3; i <= MAXA; ++i) {
    params->setSimulations(s0);
    params->setA(a[i]);
    printf("a[i]: %lf\n", a[i]);
    p0[i] = master(params);
    pMAP0[i] = *std::max_element(p0[i].begin(), p0[i].end());
  }

  vector<int> convergeGlobal(MAXA + 1, 0);
  vector<double> P;
  int bits = params->realization().realization().bitCount();
  for (int s = 1; s < (int)sims.size(); ++s) {
    int s1 = sims[s];
    printf("\ns: %d\n", s1);
    params->setSimulations(s1);
    vector<double> p1[MAXA + 1];
    vector<double> pMAP1(MAXA + 1, 0);

    for (int i = MAXA; i >= 3; --i) {
      printf("\ns: %d a: %.10lf\n", s1, a[i]);
      params->setA(a[i]);
      double converge = true;
      p1[i] = master(params);
      pMAP1[i] = *std::max_element(p1[i].begin(), p1[i].end());
      double delta = dabs(pMAP1[i] - pMAP0[i]) / pMAP1[i];
      printf("c: %lf\n", delta);
      if (delta > c) converge = false;
      int pos = 0;
      for (int j = 0; j < (int)p1[i].size(); ++j) {
        if (dabs(p1[i][j] - p0[i][j]) > c) converge = false;
        if (p1[i][j] > 0) pos++;
      }
      if (pos == 0) converge = false;
      if (converge) {
        convergeGlobal[i]++;
        printf("Converged for n=%d a=%lf\n", s0, a[i]);
        if (convergeGlobal[i] > 0) break;
      } else {
        convergeGlobal[i] = 0;
        printf("Not converged.\n");
      }
    }

    bool done = false;
    for (int i = MAXA; i >= 3; --i) {
      if (convergeGlobal[i] > 0) {
        P = p0[i];
        params->setSimulations(s0);
        params->setA(a[i]);
        done = true;
        break;
      }
    }

    s0 = s1;
    for (int i = 3; i <= MAXA; ++i) p0[i] = p1[i];
    pMAP0 = pMAP1;
    if (done) break;
    P = p1[MAXA];
  }
  return P;
}
*/
}
