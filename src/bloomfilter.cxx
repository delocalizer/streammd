#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <tuple>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <xxhash_cpp/include/xxhash.hpp>

#include "bloomfilter.h"

using namespace bloomfilter;

BloomFilter::BloomFilter(uint64_t n, float p):  n_{n}, p_{p} {
  std::tie(m_, k_) = m_k_min(n_, p_);
  bitset = std::unique_ptr<boost::dynamic_bitset<>>(new boost::dynamic_bitset<>(m_));
}

// Check if item is present.
bool BloomFilter::operator&(const std::string& item) {
  return false;
}

// Add the item; return false if it was already present otherwise true.
bool BloomFilter::operator|=(const std::string& item) {
  uint64_t hashes[k_];
  hash(item, hashes);
  for (int i = 0; i < k_; i++) {
    std::cout << hashes[i] << std::endl;
  }
  return false;
}

// Return the memory-optimal bitset size m and number of hash functions k
// for a given number of items to store n and target false positive rate p.
std::tuple<uint64_t, int> BloomFilter::m_k_min(uint64_t n, float p) {
  uint64_t m = std::ceil(n * - std::log(p) / std::pow(std::log(2), 2));
  int k = std::ceil(std::log(2) * m / n);
  return { std::tuple<uint64_t, int> { m, k } };
}

// Generate k hash values for the item.
void BloomFilter::hash(const std::string& item, uint64_t* buf) {
  uint64_t a { xxh::xxhash3<64>(item, seed1_) };
  uint64_t b { xxh::xxhash3<64>(item, seed2_) };
  for (uint64_t i = 0; i < k_; i++) {
    buf[i] = a;
    a += b;
    b += i;
  }
}
