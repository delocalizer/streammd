#ifndef STREAMMD_BLOOMFILTER_H_
#define STREAMMD_BLOOMFILTER_H_

#include <cstdint>
#include <string>
#include <tuple>

#include <dynamic_bitset/dynamic_bitset.hpp>

namespace bloomfilter {

// A simple but fast Bloom filter.
class BloomFilter {

 public:

  BloomFilter(double p, uint64_t n);
  BloomFilter(double p, uint64_t m, size_t k);

  bool contains(const std::string& item);
  bool add(const std::string& item);

  const double& p() const { return p_; }
  const uint64_t& n() const { return n_; }
  const uint64_t& m() const { return m_; }
  const size_t& k() const { return k_; }
  const bool& mpow2() const { return mpow2_; }

  uint64_t count_estimate();

  static BloomFilter fromMemSpec(double p, std::string memspec, bool mpow2=false);
  static uint64_t capacity(double p, uint64_t m, size_t k);
  static uint64_t memspec_to_bytes(std::string memspec, bool mpow2=false);
  static std::tuple<uint64_t, size_t> m_k_min(double p, uint64_t n);

 private:

  // any 2 primes "should" do
  static const uint64_t seed1_ { 43 };
  static const uint64_t seed2_ { 9967 };

  const double p_;
  
  uint64_t n_;
  uint64_t m_;
  size_t k_;
  // if m_ is a power of 2 we can use a bitmask instead of % for hash addressing
  bool mpow2_;
  uint64_t mask_;

  void initialize();
  void hash(const std::string& item, uint64_t* buf);
  std::unique_ptr<sul::dynamic_bitset<>> bitset;
  
};

}
#endif // STREAMMD_BLOOMFILTER_H_
