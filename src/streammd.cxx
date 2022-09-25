#include <cmath>
#include <iostream>
#include <string>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <xxhash_cpp/include/xxhash.hpp>
#include "bloomfilter.h"

int main(int argc, char* argv[]) {
  // initializes all to 0
  bloomfilter::BloomFilter bf(1000, 0.1);
  std::string item { "foo bar baz" };
  std::cout << "bf & item before add: " << (bf & item) << std::endl;
  bf |= item; 
  std::cout << "bf & item after add: " << (bf & item) << std::endl;
}
