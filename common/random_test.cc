#include "random.h"

#include <iostream>
#include <gtest/gtest.h>

#include <vector>

namespace {
using cmplx::common::Random;

class RandomTest : public testing::Test {};

TEST_F(RandomTest, CorrelationsTest) {
  Random r(1LL * time(NULL) * getpid());
  const int N = 1e7;
  std::vector<double> x;
  double avg = 0;
  for (int i = 0; i < N + 30; ++i) {
    x.push_back(r.prandReal01());
    if (i < N)
      avg += x[i];
  }
  avg /= N;
  std::vector<double> eps;
  for (int n = 1; n < N; ++n) {
    double cur = 0;
    for (int i = 0; i < 10; ++i) {
      cur += x[i] * x[i + n];
    }
    eps.push_back(cur / N - avg * avg);
  }

  double err = pow(N, -0.5);
  for (double e : eps) {
    if (e < 0)
      e *= -1;
    EXPECT_LT(e, err);
  }
}

TEST_F(RandomTest, MomentsTest) {
  Random r(1LL * time(NULL) * getpid());
  const int N = 1e7;
  std::vector<double> x;
  for (int i = 0; i < N; ++i)
    x.push_back(r.prandReal01());
  std::vector<double> eps;
  for (int k = 1; k <= 10; ++k) {
    double powavg = 0;
    for (int i = 0; i < N; ++i)
      powavg += pow(x[i], k);
    eps.push_back(powavg / N - 1.l / (k + 1));
  }
  double err = pow(N, -0.5);
  for (double e : eps) {
    if (e < 0)
      e *= -1;
    EXPECT_LT(e, err);
  }
}

TEST_F(RandomTest, FairDiceThrow) {
  Random r(1LL * time(NULL) * getpid());
  const int N = 1e8;
  double p = 0.5;
  int pos = 0;
  for (int i = 0; i < N; ++i) {
    if (r.eventDraw(p))
      pos++;
  }
  printf("%d %d\n", pos, N - pos);
}

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
} // namespace
