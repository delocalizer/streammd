#include <cmath>

#include <bytesize/bytesize.hh>
#include <dynamic_bitset/dynamic_bitset.hpp>
#include <xxHash/xxhash.h>
#include <spdlog/spdlog.h>

#include "bloomfilter.h"

namespace bloomfilter {

// Create a Bloom filter with false-positive rate p and n items
BloomFilter::BloomFilter(double p, uint64_t n): p_{p}, n_{n} {
  std::tie(m_, k_) = m_k_min(p_, n_);
  initialize();
}

// Create a Bloom filter with false-positive rate p, bit array size m, and k hash functions.
BloomFilter::BloomFilter(double p, uint64_t m, size_t k): p_{p}, m_{m}, k_{k} {
  n_ = capacity(p_, m_, k_);
  initialize();
}

// Check if item is present.
bool BloomFilter::contains(const std::string& item) {
  uint64_t hashes[k_];
  hash(item, hashes);
  for (size_t i = 0; i < k_; i++) {
    uint64_t pos { mpow2_ ? hashes[i] & mask_ : hashes[i] % m_ };
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
    uint64_t pos { mpow2_ ? hashes[i] & mask_ : hashes[i] % m_ };
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

// Create a Bloom filter from memory spec and false-positive rate p. For
// simplicity we use fixed k==10 and infer the capacity n. This is more useful
// in many cases than requiring n and p then inferring m and k. When mpow2 is
// true and the interpreted value of m is not a power of two, the nearest power
// of two less than the interpreted value is used.
BloomFilter BloomFilter::fromMemSpec(double p, std::string memspec, bool mpow2) {
  auto m = 8 * memspec_to_bytes(memspec, mpow2);
  return BloomFilter(p, m, 10);
}

// Return Bloom filter capacity n inferred from p, m, and k.
uint64_t BloomFilter::capacity(double p, uint64_t m, size_t k) {
  double m_dbl(m), k_dbl(k);
  return std::ceil(m_dbl / (-k_dbl / std::log(1.0 - std::exp(std::log(p) / k_dbl))));
}

// Return number of bytes parsed from memory spec. When mpow2 is true and the
// resulting number is not a power of two, return the nearest power of two less
// than the parsed value.
uint64_t BloomFilter::memspec_to_bytes(std::string memspec, bool mpow2) {
  auto m = bytesize::bytesize::parse(memspec);
  if (mpow2 && ((m & (m-1)) != 0)) {
    size_t pow = 1;
    while (pow < m) {
      pow *= 2;
    }
    m = pow >> 1;
  }
  return m;
}

// Return the memory-optimal bitset size m and number of hash functions k
// for a given number of items to store n and target false positive rate p.
std::tuple<uint64_t, size_t> BloomFilter::m_k_min(double p, uint64_t n) {
  uint64_t m = std::ceil(n * -std::log(p) / std::pow(std::log(2.0), 2.0));
  size_t k = std::ceil(std::log(2.0) * m/n);
  return { std::tuple<uint64_t, size_t> { m, k } };
}

// Initialize the bitset; set mpow2_ and mask_; log state.
void BloomFilter::initialize(){
  // dynamic_bitset constructor initializes all to 0
  bitset = std::unique_ptr<sul::dynamic_bitset<>>(new sul::dynamic_bitset(m_));
  mpow2_ = ((m_ & (m_-1)) == 0);
  mask_ = mpow2_ ? (m_-1) : 0;
  spdlog::info("BloomFilter initialized with p={0:.3g} n={1} m={2} k={3}", p_, n_, m_, k_);
  spdlog::debug("mpow2={0}; {1} will be used for hash addressing", mpow2_,
                mpow2_ ? "bit mask" : "modulus");
}

// Generate k hash values for the item.
//
// k linear combinations of just 2 independent hashes ("double hashing") has
// the same asymptotic behaviour as k independent hashes.
//
// Ref: Kirsch & Mitzenmacher (2006) https://doi.org/10.1007/11841036_42
void BloomFilter::hash(const std::string& item, uint64_t* buf) {
  auto cstr { item.c_str() };
  auto len { item.length() };
  auto a { XXH3_64bits_withSeed(cstr, len, seed1_) };
  auto b { XXH3_64bits_withSeed(cstr, len, seed2_) };
  for (size_t i = 0; i < k_; i++) {
    buf[i] = a;
    a += b;
    b += i;
  }
}

}
