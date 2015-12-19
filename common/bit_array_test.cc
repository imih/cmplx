#include <gtest/gtest.h>

class BitArrayTest : public testing::Test {
};

TEST_F(BitArrayTest, BitArraySanity) {}

int main(int argc, char* argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
