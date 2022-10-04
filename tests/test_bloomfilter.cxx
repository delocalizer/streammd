#include <tuple>

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "bloomfilter.h"

using namespace bloomfilter;

TEST_CASE("BloomFilter::m_k_min behaviour", "[static]") {
  REQUIRE(BloomFilter::m_k_min(1000000, 0.000001) == std::make_tuple(28755177, 20));
  REQUIRE(BloomFilter::m_k_min(10000000, 0.0000001) == std::make_tuple(335477051,	24));
  REQUIRE(BloomFilter::m_k_min(100000000, 0.00000001) == std::make_tuple(3834023396, 27));
  REQUIRE(BloomFilter::m_k_min(1000000000, 0.000001) == std::make_tuple(28755176136,	20));
}

TEST_CASE("BloomFilter::add missing", "[functionality]") {
  BloomFilter bf(1000, 0.001);
  auto key = "something";
  REQUIRE(bf.add(key) == true);
}

TEST_CASE("BloomFilter::add existing", "[functionality]") {
  BloomFilter bf(1000, 0.001);
  auto key = "something";
  bf.add(key);
  REQUIRE(bf.add(key) == false);
}

TEST_CASE("BloomFilter::contains missing", "[functionality]") {
  BloomFilter bf(1000, 0.001);
  auto key = "something";
  REQUIRE(bf.contains(key) == false); 
}

TEST_CASE("BloomFilter::contains existing", "[functionality]") {
  BloomFilter bf(1000, 0.001);
  auto key = "something";
  bf.add(key);
  REQUIRE(bf.contains(key) == true);
}

TEST_CASE("BloomFilter::count_estimate", "[performance]") {
  size_t n { 1000000 }, count { 0 };
  float p { 0.000001 };
  BloomFilter bf(n, p);
  for (size_t i { 0 }; i < 1000000; ++i) {
    bf.add(std::to_string(i));
    count++;
  }
  REQUIRE_THAT(float(count)/bf.count_estimate(),
               Catch::Matchers::WithinAbs(1.0, 0.001));
}

TEST_CASE("BloomFilter FNR == 0", "[performance]") {
  size_t n { 1000000 }, not_present { 0 };
  size_t imax = n;
  float p { 0.000001 };
  BloomFilter bf(n, p);
  std::vector<size_t> values(imax);
  for (size_t i { 0 }; i < imax; ++i) {
    values[i] = i;
    bf.add(std::to_string(i));
  }
  for (size_t i { 0 }; i < imax; ++i) {
    not_present += bf.contains(std::to_string(i)) ? 0 : 1;
  }
  REQUIRE(not_present == 0);
}

TEST_CASE("BloomFilter FPR bound", "[performance]") {
  size_t n { 1000000 };
  std::vector<float> ps = { 0.001, 0.0001, 0.00001, 0.000001 };
  std::vector<std::string> values(n), misses(n);
  for (size_t i { 0 }; i < n; ++i ) { values[i] = std::to_string(i); }
  for (size_t i { 0 }; i < n; ++i ) { misses[i] = std::to_string(n+i); }
  for (float p : ps) {
    BloomFilter bf(n, p);
    for (std::string value : values) {
      bf.add(value);
    }
    size_t fps { 0 };
    for (size_t i { 0 }; i < n; ++i ) { fps += bf.contains(misses[i]) ? 1 : 0; }
    auto fpr { float(fps) / n };
    // 0 <= fpr <= 2p
    REQUIRE_THAT(fpr, Catch::Matchers::WithinAbs(p, p));
  }
}
