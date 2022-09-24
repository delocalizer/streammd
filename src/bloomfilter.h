#ifndef STREAMMD_BLOOMFILTER_H_
#define STREAMMD_BLOOMFILTER_H_

#include <cstdint>
#include <string>
#include <tuple>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>

namespace bloomfilter {

class BloomFilter {

 public:

  BloomFilter(uint64_t n, float p);
  //~BloomFilter();

  bool operator&(const std::string&);
  bool operator|=(const std::string&);

  inline uint64_t n() { return n_; }
  inline float p() { return p_; }
  inline uint64_t m() { return m_; }
  inline int k() { return k_; }

  static std::tuple<uint64_t, int> m_k_min(uint64_t n, float p);

 private:

  // any 2 primes should do
  static const uint64_t seed1 = 43;
  static const uint64_t seed2 = 9967;

  const uint64_t n_;
  const float p_;

  uint64_t m_;
  int k_;

  std::unique_ptr<boost::dynamic_bitset<>> bitset;
  
};

}
#endif // STREAMMD_BLOOMFILTER_H_
