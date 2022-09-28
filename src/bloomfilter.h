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

  BloomFilter(uint64_t n, float p);

  bool contains(const std::string& item);
  bool add(const std::string& item);

  inline uint64_t n() { return n_; }
  inline float p() { return p_; }
  inline uint64_t m() { return m_; }
  inline int k() { return k_; }

  uint64_t count_estimate();

  static std::tuple<uint64_t, int> m_k_min(uint64_t n, float p);

 private:

  // any 2 primes "should" do
  static const uint64_t seed1_ { 43 };
  static const uint64_t seed2_ { 9967 };

  const uint64_t n_;
  const float p_;

  uint64_t m_;
  int k_;

  void hash(const std::string& item, uint64_t* buf);
  std::unique_ptr<sul::dynamic_bitset<>> bitset;
  
};

}
#endif // STREAMMD_BLOOMFILTER_H_
