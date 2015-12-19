#include "bit_array.h"

#include <gtest/gtest.h>

namespace {
using cmplx::common::BitArray;

class BitArrayTest : public testing::Test {};

TEST_F(BitArrayTest, BitArraySanity) {
  int n = 100;
  BitArray ba(n);
  EXPECT_EQ(n, ba.bits_num());
  for (int i = 0; i < n; ++i) {
    ba.set(i, true);
    EXPECT_EQ(i + 1, ba.bitCount());
    EXPECT_EQ(true, ba.bit(i));
  }

  for (int i = 0; i < n; ++i) {
    ba.set(i, false);
    EXPECT_EQ(false, ba.bit(i));
    EXPECT_EQ(n - i - 1, ba.bitCount());
  }
}

int main(int argc, char* argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
}  // namespace
