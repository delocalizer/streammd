#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <spdlog/spdlog.h>
#include "markdups.h"
#include "version.h"

#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

/*******************************************************************************
 * Replicate relevant tests from picard/src/test/java/picard/sam/markduplicates
 ******************************************************************************/

// (QNAME, RNAME, POS)
typedef std::tuple<std::string, std::string, size_t> SamRecordId;
typedef std::map<SamRecordId, uint16_t> SamRecordFlags;

using namespace markdups;

// Populate SAM flags map from SAM record input stream
void record_flags(std::istream& instream, SamRecordFlags& flags) {
  for (SamRecord sr; std::getline(instream, sr.buffer); ) {
    if (sr.buffer[0] != '@') {
      sr.parse();
      flags[std::make_tuple(sr.qname(), sr.rname(), sr.pos())] = sr.flag();
    }
  }
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testBulkFragmentsNoDuplicates[0]",
  "[picard tests]"){
  std::ifstream test_input, expected_output;
  SamRecordFlags marked_flags, expected_flags;
  test_input.open("resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkFragmentsNoDuplicates[0]/input.sam");
  expected_output.open("resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkFragmentsNoDuplicates[0]/output.sam");
  record_flags(expected_output, expected_flags);

  std::ostringstream test_output;
  bloomfilter::BloomFilter bf(0.000001, 1000000);
  std::vector<std::string> cli_args { "dummy", "args" };
  process_input_stream(test_input, test_output, bf, cli_args, 1);
  auto outlines { std::istringstream(test_output.str()) };
  record_flags(outlines, marked_flags);

  for (auto const& [key, val]: expected_flags) {
    INFO(fmt::format("SamRecordId: {}:{}:{}",
                     std::get<0>(key), std::get<1>(key), std::get<2>(key)));
    CHECK(marked_flags[key] == val);
  }
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testBulkFragmentsWithDuplicates[0]",
  "[picard tests]"){
  std::ifstream test_input, expected_output;
  SamRecordFlags marked_flags, expected_flags;
  test_input.open("resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkFragmentsWithDuplicates[0]/input.sam");
  expected_output.open("resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkFragmentsWithDuplicates[0]/output.sam");
  record_flags(expected_output, expected_flags);

  std::ostringstream test_output;
  bloomfilter::BloomFilter bf(0.000001, 1000000);
  std::vector<std::string> cli_args { "dummy", "args" };
  process_input_stream(test_input, test_output, bf, cli_args, 1);
  auto outlines { std::istringstream(test_output.str()) };
  record_flags(outlines, marked_flags);

  for (auto const& [key, val]: expected_flags) {
    INFO(fmt::format("SamRecordId: {}:{}:{}",
                     std::get<0>(key), std::get<1>(key), std::get<2>(key)));
    CHECK(marked_flags[key] == val);
  }
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testBulkPairsNoDuplicates[0]",
  "[picard tests]"){
  std::ifstream test_input, expected_output;
  SamRecordFlags marked_flags, expected_flags;
  test_input.open("resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkPairsNoDuplicates[0]/input.sam");
  expected_output.open("resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkPairsNoDuplicates[0]/output.sam");
  record_flags(expected_output, expected_flags);

  std::ostringstream test_output;
  bloomfilter::BloomFilter bf(0.000001, 1000000);
  std::vector<std::string> cli_args { "dummy", "args" };
  process_input_stream(test_input, test_output, bf, cli_args, 2);
  auto outlines { std::istringstream(test_output.str()) };
  record_flags(outlines, marked_flags);

  for (auto const& [key, val]: expected_flags) {
    INFO(fmt::format("SamRecordId: {}:{}:{}",
                     std::get<0>(key), std::get<1>(key), std::get<2>(key)));
    CHECK(marked_flags[key] == val);
  }
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testBulkPairsWithDuplicates[0]",
  "[picard tests]"){
  std::ifstream test_input, expected_output;
  SamRecordFlags marked_flags, expected_flags;
  test_input.open("resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkPairsWithDuplicates[0]/input.sam");
  expected_output.open("resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkPairsWithDuplicates[0]/output.sam");
  record_flags(expected_output, expected_flags);

  std::ostringstream test_output;
  bloomfilter::BloomFilter bf(0.000001, 1000000);
  std::vector<std::string> cli_args { "dummy", "args" };
  process_input_stream(test_input, test_output, bf, cli_args, 2);
  auto outlines { std::istringstream(test_output.str()) };
  record_flags(outlines, marked_flags);

  for (auto const& [key, val]: expected_flags) {
    INFO(fmt::format("SamRecordId: {}:{}:{}",
                     std::get<0>(key), std::get<1>(key), std::get<2>(key)));
    CHECK(marked_flags[key] == val);
  }
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testDifferentChromosomesInOppositeOrder",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairFirstOppositeStrandSecondUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairFirstUnmapped",
  "[picard tests][!shouldfail]"){
  INFO("streammd requires both ends of pair to match to call a duplicate");
  REQUIRE(false);
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairSecondUnmapped",
  "[picard tests][!shouldfail]"){
  INFO("streammd requires both ends of pair to match to call a duplicate");
  REQUIRE(false);
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairWithFirstEndSamePositionAndOther",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairWithSamePosition",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairWithSamePositionSameCigar",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMatePairFirstUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMatePairSecondUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testPathologicalOrderingAtTheSamePosition",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testSecondEndIsBeforeFirstInCoordinate",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testSingleMappedFragment",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testSingleMappedPair",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testSingleUnmappedFragment",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testSingleUnmappedPair",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testThreeGroupsOnDifferentChromosomesOfThreeMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testThreeMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testThreeMappedPairsWithMatchingSecondMate",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfThreeMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfTwoFragments",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfTwoMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedFragments",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePosition",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePositionDifferentStrands",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePositionDifferentStrands2",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsAndTerminalUnmappedPair",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsMatesSoftClipped",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithOppositeOrientations",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithOppositeOrientationsNumberTwo",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSoftClipping",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSoftClippingBoth",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[0]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[1]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[2]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[3]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[4]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[5]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[0]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[1]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[2]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[3]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[4]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[5]",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoUnmappedFragments",
  "[picard tests]"){
}

