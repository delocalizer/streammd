#include <cmath>
#include <cstdint>
#include <string>
#include <tuple>

#include <dynamic_bitset/dynamic_bitset.hpp>
#include <xxhash_cpp/xxhash.hpp>
#include <spdlog/spdlog.h>

#include "bloomfilter.h"

namespace bloomfilter {

// Create a Bloom filter for target maximum n items and false-positive rate p.
BloomFilter::BloomFilter(uint64_t n, float p):  n_{n}, p_{p} {
  std::tie(m_, k_) = m_k_min(n_, p_);
  // dynamic_bitset constructor initializes all to 0
  bitset = std::unique_ptr<sul::dynamic_bitset<>>(new sul::dynamic_bitset(m_));
  spdlog::info("BloomFilter initialized with n={0} p={1} m={2} k={3}", n_, p_, m_, k_);
}

// Check if item is present.
bool BloomFilter::contains(const std::string& item) {
  uint64_t hashes[k_];
  hash(item, hashes);
  for (int i = 0; i < k_; i++) {
    uint64_t pos { hashes[i] % m_ };
    if (!bitset->test(pos)) {
      return false;
    }
  }
  return true;
}

// Add the item; return false if it was already present otherwise true.
bool BloomFilter::add(const std::string& item) {
  bool added { false };
  uint64_t hashes[k_];
  hash(item, hashes);
  for (int i = 0; i < k_; i++) {
    uint64_t pos { hashes[i] % m_ };
    if (!bitset->test(pos)) {
      bitset->set(pos);
      added = true;
    }
  }
  return added;
}

// Return the estimated number of items stored.
//
// Ref: Swamidass & Baldi (2007) https://doi.org/10.1021/ci600358f
uint64_t BloomFilter::count_estimate(){
  return std::ceil((m_/k_) * - std::log(1-(float(bitset->count())/m_)));
}

// Return the memory-optimal bitset size m and number of hash functions k
// for a given number of items to store n and target false positive rate p.
std::tuple<uint64_t, int> BloomFilter::m_k_min(uint64_t n, float p) {
  uint64_t m = std::ceil(n * -std::log(p) / std::pow(std::log(2.0), 2.0));
  int k = std::ceil(std::log(2.0) * m/n);
  return { std::tuple<uint64_t, int> { m, k } };
}

// Generate k hash values for the item.
//
// k linear combinations of just 2 independent hashes ("double hashing") has
// the same asymptotic behaviour as k independent hashes.
// Ref: Kirsch & Mitzenmacher (2006) https://doi.org/10.1007/11841036_42
void BloomFilter::hash(const std::string& item, uint64_t* buf) {
  uint64_t a { xxh::xxhash3<64>(item, seed1_) };
  uint64_t b { xxh::xxhash3<64>(item, seed2_) };
  for (int i = 0; i < k_; i++) {
    buf[i] = a;
    a += b;
    b += i;
  }
}

}
