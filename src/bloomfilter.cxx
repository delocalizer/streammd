#include <cmath>

#include <bytesize/bytesize.hh>
#include <dynamic_bitset/dynamic_bitset.hpp>
#include <xxhash_cpp/xxhash.hpp>
#include <spdlog/spdlog.h>

#include "bloomfilter.h"

namespace bloomfilter {

// Create a Bloom filter with false-positive rate p and n items
BloomFilter::BloomFilter(double p, uint64_t n): p_{p}, n_{n} {
  std::tie(m_, k_) = m_k_min(p_, n_);
  // dynamic_bitset constructor initializes all to 0
  bitset = std::unique_ptr<sul::dynamic_bitset<>>(new sul::dynamic_bitset(m_));
  spdlog::info("BloomFilter initialized with p={0} n={1} m={2} k={3}", p_, n_, m_, k_);
}

// Create a Bloom filter with false-positive rate p, bit array size m, and k hash functions.
BloomFilter::BloomFilter(double p, uint64_t m, size_t k): p_{p}, m_{m}, k_{k} {
  n_ = capacity(p_, m_, k_);
  // dynamic_bitset constructor initializes all to 0
  bitset = std::unique_ptr<sul::dynamic_bitset<>>(new sul::dynamic_bitset(m_));
  spdlog::info("BloomFilter initialized with p={0} n={1} m={2} k={3}", p_, n_, m_, k_);
}

// Check if item is present.
bool BloomFilter::contains(const std::string& item) {
  uint64_t hashes[k_];
  hash(item, hashes);
  for (size_t i = 0; i < k_; i++) {
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
  for (size_t i = 0; i < k_; i++) {
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
  return std::ceil((m_/k_) * - std::log(1-(double(bitset->count())/m_)));
}

// Create a Bloom filter from memory spec (m) and false-positive rate p. For
// simplicity we use fixed k==10 and infer the capacity n. This is more useful
// in many cases than requiring n and p then inferring m and k.
BloomFilter BloomFilter::fromMemSpec(double p, std::string mem) {
  return BloomFilter(p, 8 * bytesize::bytesize::parse(mem), 10);
}

// Return Bloom filter capacity n inferred from p, m, and k
uint64_t BloomFilter::capacity(double p, uint64_t m, size_t k) {
  double m_dbl(m), k_dbl(k);
  return std::ceil(m_dbl / (-k_dbl / std::log(1.0 - std::exp(std::log(p) / k_dbl))));
}

// Return the memory-optimal bitset size m and number of hash functions k
// for a given number of items to store n and target false positive rate p.
std::tuple<uint64_t, size_t> BloomFilter::m_k_min(double p, uint64_t n) {
  uint64_t m = std::ceil(n * -std::log(p) / std::pow(std::log(2.0), 2.0));
  size_t k = std::ceil(std::log(2.0) * m/n);
  return { std::tuple<uint64_t, size_t> { m, k } };
}

// Generate k hash values for the item.
//
// k linear combinations of just 2 independent hashes ("double hashing") has
// the same asymptotic behaviour as k independent hashes.
//
// Ref: Kirsch & Mitzenmacher (2006) https://doi.org/10.1007/11841036_42
void BloomFilter::hash(const std::string& item, uint64_t* buf) {
  uint64_t a { xxh::xxhash3<64>(item, seed1_) };
  uint64_t b { xxh::xxhash3<64>(item, seed2_) };
  for (size_t i = 0; i < k_; i++) {
    buf[i] = a;
    a += b;
    b += i;
  }
}

}
