#include <cstdint>
#include <string>
#include <tuple>
#include "bloomfilter.h"

using namespace bloomfilter;

BloomFilter::BloomFilter(uint64_t n, float p):  n_(n), p_(p) {
  std::tie(m_, k_) = m_k_min(n_, p_);
}

bool BloomFilter::operator&(const std::string&) {
  return false;
}

bool BloomFilter::operator|=(const std::string&) {
  return false;
}

std::tuple<uint64_t, int> BloomFilter::m_k_min(uint64_t n, float p) {
  std::tuple<uint64_t, int> m_k { 1000, 10 };
  return m_k;
}

