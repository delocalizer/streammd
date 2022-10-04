#include "markdups.h"
#include "version.h"

#include <catch2/catch.hpp>

using namespace markdups;

TEST_CASE("markdups::SamRecord.parse", "[SamRecord]") {
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25");
  REQUIRE(sr.qname() == "HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680");
  REQUIRE(sr.flag() == 99);
  REQUIRE(sr.rname() == "chr1");
}

TEST_CASE("markdups::SamRecord.start_pos no soft clip", "[SamRecord]") {
  // this is a fwd read with POS=93578030 and CIGAR=101M
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25");
                         // POS
  REQUIRE(sr.start_pos() == 93578030);
}

TEST_CASE("markdups::SamRecord.start_pos with soft clip", "[SamRecord]") {
  // this is a fwd read with POS=725721 and CIGAR=10S32M10I34M15S
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1107:3646:4446\t99\tchr1\t725721\t30\t10S32M10I34M15S\t=\t725895\t275\tGAACACGAAAGGAATACAAGGGAATTTAATGGAATGGACTCTAATGGAATGAAATGGAATGGACTTGAATGGAATATAATGGAAGATATTAGAATGGAATA\tCCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJIJJJJIIJJJIJIJJIIJIJJJJJJEHHIJIJJJJIGHHHHHFFDEFFFEEEECDDDFEEDDDDDCDD>\tNM:i:12\tAS:i:40\tXS:i:43");
                         // POS      10S
  REQUIRE(sr.start_pos() == 725721 - 10);
}

TEST_CASE("markdups::SamRecord.end_pos with insertion and no soft clip", "[SamRecord]") {
  // this is a rev read with POS=91356686 and CIGAR=28M2I71M
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1114:14456:12207\t83\tchr1\t91356686\t60\t28M2I71M\t=\t91356467\t-318\tCTCGTATCATTTTGAACTTTGCTTTTCATTTTTTTTTTTTTAATTTATTGATTGATTGATTTATTGATCATTCTTGGGTGTTTCTCGCAGAGGGGGATTTG\t#@<A@ACCCCCDA:CCC?C?2CCCA>>DBDDD@BDFEEIIJJIJJJJJIJIIJJJIJJJJJJJJJIHGJJJJIJJJJIJJJJIIHIIIHHHHHFFFFFCCB\tNM:i:2\tAS:i:91\tXS:i:53");
                       // POS        28M  2I  71M 
  REQUIRE(sr.end_pos() == 91356686 + 28 + 0 + 71 );
}

TEST_CASE("markdups::SamRecord.end_pos with soft clip", "[SamRecord]") {
  // this is a rev read with POS=725895 and CIGAR=40M10D51M10S
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1107:3646:4446\t147\tchr1\t725895\t30\t40M10D51M10S\t=\t725721\t-275\tTGGAATGGACTCGAATGGAATGGTATGGAATGGACTCGAATGCAATGGAATGTACTCAAATGGAATGCTATGGAATTGACTCGAGTGGAATGGAATAGAAT\t;FFFFEEHFHHIIJJIJIJJIIIJIJJJJJIJJIJJIJJJJJJJJJIGIGGIGDJJGJJJIJJJJJJJJJJJJJJJJJHJJJJJJJJJHHHHHFFFFFCCC\tNM:i:17\tAS:i:40\tXS:i:35");
                       // POS      40M  10D  51M  10S
  REQUIRE(sr.end_pos() == 725895 + 40 + 10 + 51 + 10);
}

TEST_CASE("markdups::SamRecord.update_dup_status set with no previous PG", "[SamRecord]"){
}

TEST_CASE("markdups::SamRecord.update_dup_status set with previous PG last tag", "[SamRecord]"){
}

TEST_CASE("markdups::SamRecord.update_dup_status set with previous PG not last tag", "[SamRecord]"){
}

TEST_CASE("markdups::pgline with previous @PG", "[pgline]") {
  std::stringstream os;
  std::string header_prev { "@PG\tID:foo" };
  std::vector<std::string> cli_args { pgid };
  pgline(os, header_prev, cli_args);
  std::string expected {
    "@PG\tID:" + pgid + "\tPN:" + pgid + "\tCL:" + pgid + "\tVN:" +
      std::string(STREAMMD_VERSION) + "\tPP:" + "foo\n" }; 
  REQUIRE(os.str() == expected);
}

TEST_CASE("markdups::pgline with no previous @PG", "[pgline]") {
  std::stringstream os;
  std::string header_prev { "@SQ\tSN:chrMT\tLN:16569" };
  std::vector<std::string> cli_args { pgid };
  pgline(os, header_prev, cli_args);
  std::string expected {
    "@PG\tID:" + pgid + "\tPN:" + pgid + "\tCL:" + pgid + "\tVN:" +
      std::string(STREAMMD_VERSION) + "\n" }; 
  REQUIRE(os.str() == expected);
}
