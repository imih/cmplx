#include "ivector.h"
#include "ivector.cc"

#include <iostream>
#include <gtest/gtest.h>

namespace {
using cmplx::common::IVector;

class IVectorTest : public testing::Test {};

TEST_F(IVectorTest, IVectorSanity) {
  IVector<int> iv({});
  EXPECT_TRUE(iv.empty());
  EXPECT_EQ(0, iv.size());
  IVector<int> iv2({1, 2, 3});
  EXPECT_EQ(3, iv2.size());
  EXPECT_EQ(2, iv2[1]);
  iv2.set(1, -2);
  EXPECT_EQ(-2, iv2[1]);
}

int main(int argc, char* argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
}  // namespace
