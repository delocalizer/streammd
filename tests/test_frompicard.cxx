#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include "markdups.h"
#include "version.h"

#include <catch2/catch.hpp>

/*******************************************************************************
 * Replicate relevant tests from picard/src/test/java/picard/sam/markduplicates
 ******************************************************************************/
using namespace markdups;

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleUnmappedFragment",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoUnmappedFragments",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleUnmappedPair",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleMappedFragment",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedFragments",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleMappedPair",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleMappedFragmentAndSingleMappedPair",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testThreeMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSingleMappedFragmentAndTwoMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsAndTerminalUnmappedFragment",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsAndTerminalUnmappedPair",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsAndMappedSecondaryFragment",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedFragmentAndMappedPairFirstOfPairNonPrimary",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsMatesSoftClipped",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithSoftClipping",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithSoftClippingFirstOfPairOnlyNoMateCigar",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithSoftClippingBoth",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMatePairSecondUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMatePairFirstUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedFragmentAndMatePairSecondUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedFragmentAndMatePairFirstUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairAndMatePairSecondUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairAndMatePairFirstUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairAndMatePairFirstOppositeStrandSecondUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairAndMappedFragmentAndMatePairSecondUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairAndMappedFragmentAndMatePairFirstUnmapped",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithOppositeOrientations",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithOppositeOrientationsNumberTwo",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testThreeMappedPairsWithMatchingSecondMate",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairWithSamePosition",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairWithSamePositionSameCigar",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairWithSamePosition",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairWithSamePositionDifferentStrands",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairWithSamePositionDifferentStrands2",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testMappedPairWithFirstEndSamePositionAndOther",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoGroupsOnDifferentChromosomesOfTwoFragments",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoGroupsOnDifferentChromosomesOfTwoMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoGroupsOnDifferentChromosomesOfThreeMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testThreeGroupsOnDifferentChromosomesOfThreeMappedPairs",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testBulkFragmentsNoDuplicates",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testBulkPairsNoDuplicates",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testBulkFragmentsWithDuplicates",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testBulkPairsWithDuplicates",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testSecondEndIsBeforeFirstInCoordinate",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testPathologicalOrderingAtTheSamePosition",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testDifferentChromosomesInOppositeOrder",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithSupplementaryReads",
  "[picard tests]"){
}

TEST_CASE(
  "AbstractMarkDuplicatesCommandLineProgramTest.testTwoMappedPairsWithSupplementaryReadsAfterCanonical",
  "[picard tests]"){
}

TEST_CASE(
  "AsIsMarkDuplicatesTest.testSameUnclipped5PrimeOppositeStrand",
  "[picard tests]"){
}

TEST_CASE(
  "MarkDuplicatesTest.testTwoMappedPairsWithSoftClippingFirstOfPairOnly",
  "[picard tests]"){
}

