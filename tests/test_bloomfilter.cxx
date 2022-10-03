#ifndef CATCH_CONFIG_MAIN
#define CATCH_CONFIG_MAIN
#endif

#include <tuple>

#include <catch2/catch.hpp>

#include "bloomfilter.h"
#include "test_main.h"

using namespace bloomfilter;

TEST_CASE(" m_k_min ", "[bloomfilter]") {
  REQUIRE(BloomFilter::m_k_min(1000, 0.001) == std::make_tuple(14378,	10));
}
