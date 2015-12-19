#include "bit_array.h"

#include <bitset>
#include <chrono>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

using cmplx::common::BitArray;
using std::bitset;
using std::vector;
using std::chrono::system_clock;
using std::chrono::duration_cast;

typedef std::chrono::microseconds TimeT;

int main() {
  const int n = 1000;
  const int iterations = 1000;
  {  // constructing BitArray
    vector<BitArray> ba;
    auto start = system_clock::now();
    for (int i = 0; i < iterations; ++i) {
      BitArray ba(n);
    }
    std::cout << "Contructing BitArray (x" << iterations << ", size " << n
              << "): "
              << duration_cast<TimeT>(system_clock::now() - start).count()
              << " mis" << std::endl;
  }

  {  // contructing bitset
    auto start = system_clock::now();
    for (int i = 0; i < iterations; ++i) {
      bitset<n> bs;
    }
    std::cout << "Contructing bitset (x" << iterations << ", size " << n
              << "): "
              << duration_cast<TimeT>(system_clock::now() - start).count()
              << " mis" << std::endl;
  }

  BitArray ba(n);
  std::cout << "BitArray size: " << sizeof ba << std::endl;
  {  // setting bits, BitArray
    auto start = system_clock::now();
    for (int i = 0; i < n; i += 2) {
      ba.set(i, true);
    }
    std::cout << "Seting bits, BitArray (x" << n / 2 << "): "
              << duration_cast<TimeT>(system_clock::now() - start).count()
              << " mis" << std::endl;
  }

  {  // Clearing bits, BitArray
    auto start = system_clock::now();
    for (int i = 0; i < n; i += 2) {
      ba.set(i, false);
    }
    std::cout << "Clearing bits, BitArray (x" << n / 2 << "): "
              << duration_cast<TimeT>(system_clock::now() - start).count()
              << " mis" << std::endl;
  }

  bitset<n> bs;
  std::cout << "bitset size: " << sizeof bs << std::endl;
  {  // Setting bits, bitset
    auto start = system_clock::now();
    for (int i = 0; i < n; i += 2) {
      if (!bs.test(i)) bs.flip(i);
    }
    std::cout << "Seting bits, bitset (x" << n / 2 << "): "
              << duration_cast<TimeT>(system_clock::now() - start).count()
              << " mis" << std::endl;
  }

  {  // Clearing bits, bitset
    auto start = system_clock::now();
    for (int i = 0; i < n; i += 2) {
      if (bs.test(i)) bs.flip(i);
    }
    std::cout << "Clearing bits, bitset (x" << n / 2 << "): "
              << duration_cast<TimeT>(system_clock::now() - start).count()
              << " mis" << std::endl;
  }

  srand(time(NULL));
  {  // Random seting, BitArray
    auto start = system_clock::now();
    for (int i = 0; i < iterations; i += 2) {
      int r = rand() % n;
      ba.set(r, true);
    }
    std::cout << "Random seting bits, BitArray (x" << iterations / 2 << "): "
              << duration_cast<TimeT>(system_clock::now() - start).count()
              << " mis" << std::endl;
  }

  {  // Random setting, bitset
    auto start = system_clock::now();
    for (int i = 0; i < iterations; i += 2) {
      int r = rand() % n;
      if (!bs.test(r)) bs.flip(r);
    }
    std::cout << "Random seting bits, bitset (x" << iterations / 2 << "): "
              << duration_cast<TimeT>(system_clock::now() - start).count()
              << " mis" << std::endl;
  }

  int bc;
  {  // bitcount, BitArray
    auto start = system_clock::now();
    for (int i = 0; i < iterations; ++i) bc = ba.bitCount();
    std::cout << "Bitcount, bitarray (x" << iterations << "): "
              << duration_cast<TimeT>(system_clock::now() - start).count()
              << " mis" << std::endl;
  }

  {  // bitcount, bitset
    auto start = system_clock::now();
    for (int i = 0; i < iterations; ++i) bc = bs.count();
    std::cout << "Bitcount, bitset (x" << iterations << "): "
              << duration_cast<TimeT>(system_clock::now() - start).count()
              << " mis" << std::endl;
  }

  return 0;
}
