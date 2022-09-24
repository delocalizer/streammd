#ifndef STREAMMD_BLOOMFILTER_H_
#define STREAMMD_BLOOMFILTER_H_

#include <cstdint>
#include <tuple>

namespace bloomfilter {

class BloomFilter {

 public:

  BloomFilter(uint64_t n, float p);
  bool operator&(const std::string&);
  bool operator|=(const std::string&);

 private:

  // any 2 primes should do
  static const uint64_t seed1 = 43;
  static const uint64_t seed2 = 9967;

  const uint64_t n_;
  const float p_;

  uint64_t m_;
  int k_;
  
};

}
#endif // STREAMMD_BLOOMFILTER_H_
