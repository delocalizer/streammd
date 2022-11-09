#include <catch2/catch.hpp>

#include "test_util.h"


TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testBulkFragmentsNoDuplicates[0]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkFragmentsNoDuplicates[0]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkFragmentsNoDuplicates[0]/output.sam",
    1);
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testBulkFragmentsWithDuplicates[0]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkFragmentsWithDuplicates[0]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkFragmentsWithDuplicates[0]/output.sam",
    1);
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testBulkPairsNoDuplicates[0]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkPairsNoDuplicates[0]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkPairsNoDuplicates[0]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testBulkPairsWithDuplicates[0]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkPairsWithDuplicates[0]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testBulkPairsWithDuplicates[0]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testDifferentChromosomesInOppositeOrder",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testDifferentChromosomesInOppositeOrder/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testDifferentChromosomesInOppositeOrder/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairFirstOppositeStrandSecondUnmapped",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairFirstOppositeStrandSecondUnmapped/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairFirstOppositeStrandSecondUnmapped/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairFirstUnmapped",
  "[picard tests][!shouldfail]"){
  INFO("streammd requires both ends of pair to match to call a duplicate");
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairFirstUnmapped/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairFirstUnmapped/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairSecondUnmapped",
  "[picard tests][!shouldfail]"){
  INFO("streammd requires both ends of pair to match to call a duplicate");
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairSecondUnmapped/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairAndMatePairSecondUnmapped/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairWithFirstEndSamePositionAndOther",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairWithFirstEndSamePositionAndOther/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairWithFirstEndSamePositionAndOther/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairWithSamePosition",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairWithSamePosition/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairWithSamePosition/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMappedPairWithSamePositionSameCigar",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairWithSamePositionSameCigar/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMappedPairWithSamePositionSameCigar/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMatePairFirstUnmapped",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMatePairFirstUnmapped/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMatePairFirstUnmapped/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testMatePairSecondUnmapped",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMatePairSecondUnmapped/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testMatePairSecondUnmapped/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testPathologicalOrderingAtTheSamePosition",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testPathologicalOrderingAtTheSamePosition/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testPathologicalOrderingAtTheSamePosition/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testSecondEndIsBeforeFirstInCoordinate",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testSecondEndIsBeforeFirstInCoordinate/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testSecondEndIsBeforeFirstInCoordinate/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testSingleMappedFragment",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testSingleMappedFragment/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testSingleMappedFragment/output.sam",
    1);
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testSingleMappedPair",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testSingleMappedPair/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testSingleMappedPair/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testSingleUnmappedFragment",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testSingleUnmappedFragment/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testSingleUnmappedFragment/output.sam",
    1);
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testSingleUnmappedPair",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testSingleUnmappedPair/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testSingleUnmappedPair/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testThreeGroupsOnDifferentChromosomesOfThreeMappedPairs",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testThreeGroupsOnDifferentChromosomesOfThreeMappedPairs/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testThreeGroupsOnDifferentChromosomesOfThreeMappedPairs/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testThreeMappedPairs",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testThreeMappedPairs/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testThreeMappedPairs/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testThreeMappedPairsWithMatchingSecondMate",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testThreeMappedPairsWithMatchingSecondMate/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testThreeMappedPairsWithMatchingSecondMate/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfThreeMappedPairs",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfThreeMappedPairs/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfThreeMappedPairs/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfTwoFragments",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfTwoFragments/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfTwoFragments/output.sam",
    1);
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfTwoMappedPairs",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfTwoMappedPairs/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoGroupsOnDifferentChromosomesOfTwoMappedPairs/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedFragments",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedFragments/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedFragments/output.sam",
    1);
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePosition",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePosition/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePosition/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePositionDifferentStrands",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePositionDifferentStrands/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePositionDifferentStrands/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePositionDifferentStrands2",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePositionDifferentStrands2/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairWithSamePositionDifferentStrands2/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairs",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairs/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairs/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsAndTerminalUnmappedPair",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsAndTerminalUnmappedPair/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsAndTerminalUnmappedPair/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsMatesSoftClipped",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsMatesSoftClipped/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsMatesSoftClipped/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithOppositeOrientations",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithOppositeOrientations/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithOppositeOrientations/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithOppositeOrientationsNumberTwo",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithOppositeOrientationsNumberTwo/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithOppositeOrientationsNumberTwo/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSoftClipping",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSoftClipping/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSoftClipping/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSoftClippingBoth",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSoftClippingBoth/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSoftClippingBoth/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[0]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[0]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[0]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[1]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[1]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[1]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[2]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[2]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[2]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[3]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[3]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[3]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[4]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[4]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[4]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[5]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[5]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReadsAfterCanonical[5]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[0]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[0]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[0]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[1]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[1]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[1]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[2]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[2]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[2]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[3]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[3]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[3]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[4]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[4]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[4]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[5]",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[5]/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoMappedPairsWithSupplementaryReads[5]/output.sam");
}

TEST_CASE(
  "MarkDuplicatesTestQueryNameSorted.testTwoUnmappedFragments",
  "[picard tests]"){
  test_streammd(
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoUnmappedFragments/input.sam",
    "resources/picard_tests/MarkDuplicatesTestQueryNameSorted.testTwoUnmappedFragments/output.sam",
    1);
}

