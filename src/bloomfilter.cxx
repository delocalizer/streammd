#include <cstdint>
#include <string>
#include <tuple>
#include "bloomfilter.h"

using namespace bloomfilter;

BloomFilter::BloomFilter(uint64_t n, float p):  n_(n), p_(p) {
}

bool BloomFilter::operator&(const std::string&) {
  return false;
}

bool BloomFilter::operator|=(const std::string&) {
  return false;
}
