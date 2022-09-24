#include <iostream>
#include <cmath>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <xxhash_cpp/include/xxhash.hpp>
#include "bloomfilter.h"

int main(int argc, char* argv[]) {
  // initializes all to 0
  bloomfilter::BloomFilter bf(1000, 0.1);
  std::cout << "bf.m: " << bf.m() << std::endl;
  boost::dynamic_bitset<> b1(pow(1024, 2));
  std::cout << "test(0) before set: " << b1.test(0) << std::endl;
  b1.set(0);
  std::cout << "test(0) after set: " << b1.test(0) << std::endl;
  std::array<int, 4> input {322, 2137, 42069, 65536};
  std::string input2 = "foo bar baz";
  uint64_t seed1 = 41;
  uint64_t seed2 = 43;
  xxh::hash_t<64> hash1 = xxh::xxhash3<64>(input2, seed1);
  xxh::hash_t<64> hash2 = xxh::xxhash3<64>(input2, seed2);
  std::cout << "hash1 of input is " << hash1 << std::endl;
  std::cout << "hash1 of input is " << hash2 << std::endl;
}
