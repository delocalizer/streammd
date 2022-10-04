#include "markdups.h"
#include "version.h"

#include <catch2/catch.hpp>

using namespace markdups;

TEST_CASE("markdups::pgline with previous @PG") {
  std::stringstream os;
  std::string header_prev { "@PG\tID:foo" };
  std::vector<std::string> cli_args { pgid };
  pgline(os, header_prev, cli_args);
  std::string expected {
    "@PG\tID:" + pgid + "\tPN:" + pgid + "\tCL:" + pgid + "\tVN:" +
      std::string(STREAMMD_VERSION) + "\tPP:" + "foo\n" }; 
  REQUIRE(os.str() == expected);
}
