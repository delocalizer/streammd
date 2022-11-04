#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include "markdups.h"
#include "version.h"

#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

/*******************************************************************************
 * Replicate relevant tests from picard/src/test/java/picard/sam/markduplicates
 ******************************************************************************/

using namespace markdups;

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleUnmappedFragment",
  "[picard tests]"){
  std::stringstream in {
    "READ1\t4\t*\t0\t0\t*\t*\t0\t10\tAAAAAAAAAA\t++++++++++\n"
  };
  std::stringstream out;
  bloomfilter::BloomFilter bf(0.001, 1000);
  std::vector<std::string> cli_args { "dummy", "args" };
  auto result = process_input_stream(in, out, bf, cli_args, 1);
  std::string line;
  std::getline(out, line); // consume @PG line
  std::getline(out, line);
  SamRecord sr1(line);
  CHECK(!(sr1.flag() & flag_duplicate));
  CHECK(result.templates == 1);
  CHECK(result.templates_marked_duplicate == 0);
  CHECK(result.alignments == 1);
  CHECK(result.alignments_marked_duplicate == 0);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoUnmappedFragments",
  "[picard tests]"){
  std::stringstream in {
    "READ1\t4\t*\t0\t0\t*\t*\t0\t10\tAAAAAAAAAA\t++++++++++\n"
    "READ2\t4\t*\t0\t0\t*\t*\t0\t10\tAAAAAAAAAA\t++++++++++\n"
  };
  std::stringstream out;
  bloomfilter::BloomFilter bf(0.001, 1000);
  std::vector<std::string> cli_args { "dummy", "args" };
  auto result = process_input_stream(in, out, bf, cli_args, 1);
  std::string line;
  std::getline(out, line); // consume @PG line
  std::getline(out, line);
  SamRecord sr1(line);
  std::getline(out, line);
  SamRecord sr2(line);
  CHECK(!(sr1.flag() & flag_duplicate));
  CHECK(!(sr2.flag() & flag_duplicate));
  CHECK(result.templates == 2);
  CHECK(result.templates_marked_duplicate == 0);
  CHECK(result.alignments == 2);
  CHECK(result.alignments_marked_duplicate == 0);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleUnmappedPair",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleMappedFragment",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedFragments",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleMappedPair",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleMappedFragmentAndSingleMappedPair",
  "[picard tests][!shouldfail]"){
  INFO("streammd checks either one read or two for all templates, not both.");
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairs",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testThreeMappedPairs",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleMappedFragmentAndTwoMappedPairs",
  "[picard tests][!shouldfail]"){
  INFO("streammd checks either one read or two for all templates, not both.");
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsAndTerminalUnmappedFragment",
  "[picard tests][!shouldfail]"){
  INFO("streammd checks either one read or two for all templates, not both.");
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsAndTerminalUnmappedPair",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testOpticalDuplicateFinding",
  "[picard tests][!shouldfail]"){
  INFO("streammd cannot distinguish optical dupes from library dupes.")
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testOpticalDuplicateClusterSamePositionNoOpticalDuplicates",
  "[picard tests][!shouldfail]"){
  INFO("streammd cannot distinguish optical dupes from library dupes.")
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testOpticalDuplicateClusterSamePositionNoOpticalDuplicatesWithinPixelDistance",
  "[picard tests][!shouldfail]"){
  INFO("streammd cannot distinguish optical dupes from library dupes.")
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testOpticalDuplicateClusterSamePositionOneOpticalDuplicatesWithinPixelDistance",
  "[picard tests][!shouldfail]"){
  INFO("streammd cannot distinguish optical dupes from library dupes.")
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testOpticalDuplicateClusterOneEndSamePositionOneCluster",
  "[picard tests][!shouldfail]"){
  INFO("streammd cannot distinguish optical dupes from library dupes.")
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testManyOpticalDuplicateClusterOneEndSamePositionOneCluster",
  "[picard tests][!shouldfail]"){
  INFO("streammd cannot distinguish optical dupes from library dupes.")
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsAndMappedSecondaryFragment",
  "[picard tests][!shouldfail]"){
  INFO("streammd checks either one read or two for all templates, not both.");
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedFragmentAndMappedPairFirstOfPairNonPrimary",
  "[picard tests][!shouldfail]"){
  INFO("streammd checks either one read or two for all templates, not both.");
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsMatesSoftClipped",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithSoftClipping",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithSoftClippingFirstOfPairOnlyNoMateCigar",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithSoftClippingBoth",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMatePairSecondUnmapped",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMatePairFirstUnmapped",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedFragmentAndMatePairSecondUnmapped",
  "[picard tests][!shouldfail]"){
  INFO("streammd checks either one read or two for all templates, not both.");
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedFragmentAndMatePairFirstUnmapped",
  "[picard tests][!shouldfail]"){
  INFO("streammd checks either one read or two for all templates, not both.");
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairAndMatePairSecondUnmapped",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairAndMatePairFirstUnmapped",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairAndMatePairFirstOppositeStrandSecondUnmapped",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairAndMappedFragmentAndMatePairSecondUnmapped",
  "[picard tests][!shouldfail]"){
  INFO("streammd checks either one read or two for all templates, not both.");
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairAndMappedFragmentAndMatePairFirstUnmapped",
  "[picard tests][!shouldfail]"){
  INFO("streammd checks either one read or two for all templates, not both.");
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithOppositeOrientations",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithOppositeOrientationsNumberTwo",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testThreeMappedPairsWithMatchingSecondMate",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairWithSamePosition",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairWithSamePositionSameCigar",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairWithSamePosition",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairWithSamePositionDifferentStrands",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairWithSamePositionDifferentStrands2",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairWithFirstEndSamePositionAndOther",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoGroupsOnDifferentChromosomesOfTwoFragments",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoGroupsOnDifferentChromosomesOfTwoMappedPairs",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoGroupsOnDifferentChromosomesOfThreeMappedPairs",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testThreeGroupsOnDifferentChromosomesOfThreeMappedPairs",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testBulkFragmentsNoDuplicates",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testBulkPairsNoDuplicates",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testBulkFragmentsWithDuplicates",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testBulkPairsWithDuplicates",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testStackOverFlowPairSetSwap",
  "[picard tests]"){
  SUCCEED("Not an issue in streammd");
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSecondEndIsBeforeFirstInCoordinate",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testPathologicalOrderingAtTheSamePosition",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testDifferentChromosomesInOppositeOrder",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testOpticalDuplicateClustersAddingSecondEndFirstSameCoordinate",
  "[picard tests][!shouldfail]"){
  INFO("streammd cannot distinguish optical dupes from library dupes.")
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testOpticalDuplicateBehaviorNullRegex",
  "[picard tests][!shouldfail]"){
  INFO("streammd cannot distinguish optical dupes from library dupes.")
  REQUIRE(false);
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithSupplementaryReads",
  "[picard tests]"){
  // TODO
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithSupplementaryReadsAfterCanonical",
  "[picard tests]"){
  // TODO
}

TEST_CASE(
  "AsIsMarkDuplicatesTest.testSameUnclipped5PrimeOppositeStrand",
  "[picard tests]"){
  // TODO
}

TEST_CASE(
  "AsIsMarkDuplicatesTest.testQueryGroupedInput",
  "[picard tests]"){
  SUCCEED("QNAME grouped input is the only kind that streammd handles anyway");
}

TEST_CASE(
  "MarkDuplicatesTest.testTwoMappedPairsWithSoftClippingFirstOfPairOnly",
  "[picard tests]"){
  //TODO
}

TEST_CASE(
  "MarkDuplicatesTest.pgRecordChainingTest",
  "[picard tests]"){
  SUCCEED("PG record operations are already tested in test_markdups");
}

TEST_CASE(
  "MarkDuplicatesTest.testOpticalDuplicateDetection",
  "[picard tests][!shouldfail]"){
  INFO("streammd cannot distinguish optical dupes from library dupes.")
  REQUIRE(false);
}

TEST_CASE(
  "MarkDuplicatesTest.testWithBarcodeFragmentDuplicate",
  "[picard tests][!shouldfail]"){
  INFO("streammd currently has no specific handling of barcodes.")
  REQUIRE(false);
}

TEST_CASE(
  "MarkDuplicatesTest.testWithBarcodeDuplicate",
  "[picard tests][!shouldfail]"){
  INFO("streammd currently has no specific handling of barcodes.")
  REQUIRE(false);
}

TEST_CASE(
  "MarkDuplicatesTest.testWithBarcodeComplex",
  "[picard tests][!shouldfail]"){
  INFO("streammd currently has no specific handling of barcodes.")
  REQUIRE(false);
}

TEST_CASE(
  "MarkDuplicatesTest.testWithIndividualReadBarcodes",
  "[picard tests][!shouldfail]"){
  INFO("streammd currently has no specific handling of barcodes.")
  REQUIRE(false);
}

TEST_CASE(
  "MarkDuplicatesTest.testWithUMIs",
  "[picard tests][!shouldfail]"){
  INFO("streammd currently has no specific handling of UMIs.")
  REQUIRE(false);
}
