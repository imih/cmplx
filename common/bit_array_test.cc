#include "./bit_array.h"

#include <gtest/gtest.h>
#include <vector>

namespace {
using cmplx::common::BitArray;
class BitArrayTest : public testing::Test {};

TEST_F(BitArrayTest, BitArraySanity) {
  int n = 100;
  BitArray ba(n);
  EXPECT_EQ(n, ba.bits_num());
  std::vector<int> pos;
  for (int i = 0; i < n; ++i) {
    ba.set(i, true);
    EXPECT_EQ(i + 1, ba.bitCount());
    EXPECT_EQ(true, ba.bit(i));
    pos.push_back(i);
    EXPECT_EQ(pos, ba.positions());
  }

  for (int i = 0; i < n; ++i) {
    ba.set(i, false);
    EXPECT_EQ(false, ba.bit(i));
    EXPECT_EQ(n - i - 1, ba.bitCount());
  }
}

TEST_F(BitArrayTest, BitArrayOperations) {
  BitArray zero = BitArray::zeros(10);
  BitArray one = BitArray::ones(10);

  EXPECT_EQ(10, zero.bits_num());

  EXPECT_EQ(one, zero | one);
  EXPECT_EQ(one, one | zero);
  EXPECT_EQ(one, one | one);
  EXPECT_EQ(zero, zero | zero);

  EXPECT_EQ(zero, zero & one);
  EXPECT_EQ(zero, one & zero);
  EXPECT_EQ(one, one & one);
  EXPECT_EQ(zero, zero & zero);

  // XNOR
  EXPECT_EQ(zero, zero % one);
  EXPECT_EQ(zero, one % zero);
  EXPECT_EQ(one % one, one);
  EXPECT_EQ(zero % zero, one);

  // ==
  EXPECT_EQ(true, zero == zero);
  EXPECT_EQ(true, one == one);
  EXPECT_EQ(false, zero == one);
  EXPECT_EQ(false, BitArray::ones(1) == one);
}

// TODO test XnorSimilarity
// TODO test JaccardSimilarity

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
}  // namespace
