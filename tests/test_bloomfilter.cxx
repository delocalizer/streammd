#ifndef CATCH_CONFIG_MAIN
#define CATCH_CONFIG_MAIN
#endif

#include <tuple>

#include <catch2/catch.hpp>

#include "bloomfilter.h"

using namespace bloomfilter;

TEST_CASE("m_k_min", "[static]") {
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
