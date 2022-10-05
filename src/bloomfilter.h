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

  BloomFilter(uint64_t n, double p);
  BloomFilter(uint64_t m, double p, size_t k);

  bool contains(const std::string& item);
  bool add(const std::string& item);

  const double& p() const { return p_; }
  const size_t& k() const { return k_; }
  const uint64_t& m() const { return m_; }
  const uint64_t& n() const { return n_; }

  uint64_t count_estimate();

  static BloomFilter fromMemSpec(std::string memspec, double p);
  static uint64_t capacity(uint64_t m, size_t k, double p);
  static std::tuple<uint64_t, size_t> m_k_min(uint64_t n, double p);

 private:

  // any 2 primes "should" do
  static const uint64_t seed1_ { 43 };
  static const uint64_t seed2_ { 9967 };

  const double p_;

  size_t k_;
  uint64_t m_;
  uint64_t n_;

  void hash(const std::string& item, uint64_t* buf);
  std::unique_ptr<sul::dynamic_bitset<>> bitset;
  
};

}
#endif // STREAMMD_BLOOMFILTER_H_
