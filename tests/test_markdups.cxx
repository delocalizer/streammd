#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include "markdups.h"
#include "version.h"

#include <catch2/catch.hpp>

using namespace markdups;

// helper for testing, reimplements lines from process_qname_group
std::string ends_str(
    std::deque<std::string> ends,
    size_t reads_per_template) {
  return (reads_per_template == 1)
    ? ends.front()
    : ends.front() + '_' + ends.back();
}

TEST_CASE("markdups::SamRecord.parse", "[SamRecord]") {
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25");
  CHECK(sr.qname() == "HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680");
  CHECK(sr.flag() == 99);
  CHECK(sr.rname() == "chr1");
}

TEST_CASE("markdups::SamRecord.start_pos no soft clip", "[SamRecord]") {
  // this is a fwd read with POS=93578030 and CIGAR=101M
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25");
                         // POS
  CHECK(sr.start_pos() == 93578030);
}

TEST_CASE("markdups::SamRecord.start_pos with soft clip", "[SamRecord]") {
  // this is a fwd read with POS=725721 and CIGAR=10S32M10I34M15S
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1107:3646:4446\t99\tchr1\t725721\t30\t10S32M10I34M15S\t=\t725895\t275\tGAACACGAAAGGAATACAAGGGAATTTAATGGAATGGACTCTAATGGAATGAAATGGAATGGACTTGAATGGAATATAATGGAAGATATTAGAATGGAATA\tCCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJIJJJJIIJJJIJIJJIIJIJJJJJJEHHIJIJJJJIGHHHHHFFDEFFFEEEECDDDFEEDDDDDCDD>\tNM:i:12\tAS:i:40\tXS:i:43");
                         // POS      10S
  CHECK(sr.start_pos() == 725721 - 10);
}

TEST_CASE("markdups::SamRecord.end_pos with insertion and no soft clip", "[SamRecord]") {
  // this is a rev read with POS=91356686 and CIGAR=28M2I71M
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1114:14456:12207\t83\tchr1\t91356686\t60\t28M2I71M\t=\t91356467\t-318\tCTCGTATCATTTTGAACTTTGCTTTTCATTTTTTTTTTTTTAATTTATTGATTGATTGATTTATTGATCATTCTTGGGTGTTTCTCGCAGAGGGGGATTTG\t#@<A@ACCCCCDA:CCC?C?2CCCA>>DBDDD@BDFEEIIJJIJJJJJIJIIJJJIJJJJJJJJJIHGJJJJIJJJJIJJJJIIHIIIHHHHHFFFFFCCB\tNM:i:2\tAS:i:91\tXS:i:53");
                       // POS        28M  2I  71M 
  CHECK(sr.end_pos() == 91356686 + 28 + 0 + 71 );
}

TEST_CASE("markdups::SamRecord.end_pos with soft clip", "[SamRecord]") {
  // this is a rev read with POS=725895 and CIGAR=40M10D51M10S
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1107:3646:4446\t147\tchr1\t725895\t30\t40M10D51M10S\t=\t725721\t-275\tTGGAATGGACTCGAATGGAATGGTATGGAATGGACTCGAATGCAATGGAATGTACTCAAATGGAATGCTATGGAATTGACTCGAGTGGAATGGAATAGAAT\t;FFFFEEHFHHIIJJIJIJJIIIJIJJJJJIJJIJJIJJJJJJJJJIGIGGIGDJJGJJJIJJJJJJJJJJJJJJJJJHJJJJJJJJJHHHHHFFFFFCCC\tNM:i:17\tAS:i:40\tXS:i:35");
                       // POS      40M  10D  51M  10S
  CHECK(sr.end_pos() == 725895 + 40 + 10 + 51 + 10);
}

TEST_CASE("markdups::SamRecord.update_dup_status set with no previous PG", "[SamRecord]"){
  // FLAG=99
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25");
  sr.update_dup_status(true);
  CHECK(sr.flag() == 1123);
  CHECK(sr.buffer == "HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t1123\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25\tPG:Z:streammd\n");
}

TEST_CASE("markdups::SamRecord.update_dup_status set with previous PG last tag", "[SamRecord]"){
  // FLAG=99
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25\tPG:Z:foo");
  sr.update_dup_status(true);
  CHECK(sr.flag() == 1123);
  CHECK(sr.buffer == "HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t1123\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25\tPG:Z:streammd\n");
}

TEST_CASE("markdups::SamRecord.update_dup_status set with previous PG not last tag", "[SamRecord]"){
  // FLAG=99
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tPG:Z:foo\tXS:i:25");
  sr.update_dup_status(true);
  CHECK(sr.flag() == 1123);
  CHECK(sr.buffer == "HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t1123\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tPG:Z:streammd\tXS:i:25\n");
}

TEST_CASE("markdups::SamRecord.update_dup_status unset with no previous PG", "[SamRecord]"){
  // FLAG=1123
  SamRecord sr("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t1123\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25");
  sr.update_dup_status(false);
  CHECK(sr.flag() == 99);
  CHECK(sr.buffer == "HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25\tPG:Z:streammd\n");
}

TEST_CASE("markdups::pgline with previous @PG", "[pgline]") {
  std::stringstream os;
  std::string header_prev { "@PG\tID:foo" };
  std::vector<std::string> cli_args { pgid };
  pgline(os, header_prev, cli_args);
  std::string expected {
    "@PG\tID:" + pgid + "\tPN:" + pgid + "\tCL:" + pgid + "\tVN:" +
      std::string(STREAMMD_VERSION) + "\tPP:" + "foo\n" }; 
  CHECK(os.str() == expected);
}

TEST_CASE("markdups::pgline with no previous @PG", "[pgline]") {
  std::stringstream os;
  std::string header_prev { "@SQ\tSN:chrMT\tLN:16569" };
  std::vector<std::string> cli_args { pgid };
  pgline(os, header_prev, cli_args);
  std::string expected {
    "@PG\tID:" + pgid + "\tPN:" + pgid + "\tCL:" + pgid + "\tVN:" +
      std::string(STREAMMD_VERSION) + "\n" }; 
  CHECK(os.str() == expected);
}

TEST_CASE("markdups::template_ends pair", "[template_ends]"){
  // Confirm that ends of a pair are calculated as expected
  SamRecord pair1_r1("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25");
  SamRecord pair1_r2("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t147\tchr1\t93578228\t60\t101M\t=\t93578030\t-299\tAACAACAACAAAAAATTTGGTATTTCTAAGATGAAATGGCCAAGGCTTTCTAGTCAATTGGATTTAGAGTAAAGGAGACTATAGAAGATTACTAAGCTATA\tBDDDDDDEDDFHHHHHHJIJIJJJJJIJJJIJJJJJJJJJIIJJIJJJIIIJJIJJJJJIJJJJIJJJJJJJJJJJJJJJIJJJJJJJHHHHHFFFFFCCB\tNM:i:0\tAS:i:101\tXS:i:21");
  std::vector<SamRecord> qn1 { pair1_r1, pair1_r2 };
  auto ends1 { ends_str(template_ends(qn1), 2) };
  CHECK(ends1 == "chr1F93578030_chr1R93578329");
}

TEST_CASE("markdups::template_ends pair with one end unmapped", "[template_ends]"){
  // Confirm that ends of pair with first end unmapped are calculated as expected
  SamRecord pair1_r1("HWI-ST1213:151:C1DTBACXX:2:1104:12719:81398\t69\tchr1\t66264\t0\t*\t=\t66264\t0\tTATATGATAATATATTATTCTATAATATATTATAATTATATTATTATATTTTATTATAGGTTATATATATTATAATTATATTATATTATTATAATATGTAT\tCCCFFFFFHHHHHJJJIJJJJJJJJJJJJJJJJJJJJJJIJJJJJJJJJJJJIIJJJJIJHHJGJJJJJJJJJJJHIJJJJJJJJJJJJJJJJJIJIIJJH\tPG:Z:MarkDuplicates\tAS:i:0\tXS:i:0");
  SamRecord pair1_r2("HWI-ST1213:151:C1DTBACXX:2:1104:12719:81398\t1161\tchr1\t66264\t0\t52S10M1D39M\t=\t66264\t0\tAATTCTGTACAATAATATAATATAATTATAATATGTAATATATGATTATAATAATATAATATAATTATAATATATAATATATAATATAAATATAATATAAT\t@BCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJIJIJJJHIJJJJJJIIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIIGJJJJJJJIJJJH\tSA:Z:chrX,77138559,+,10S66M25S,0,7;\tMD:Z:10^A7T30A0\tPG:Z:MarkDuplicates\tNM:i:3\tAS:i:36\tXS:i:35");
  std::vector<SamRecord> qn1 { pair1_r1, pair1_r2 };
  auto ends1 { ends_str(template_ends(qn1), 2) };
  // First read is unmapped but 'unmapped' sorts to last
  CHECK(ends1 == "chr1F66212_" + unmapped);
}

TEST_CASE("markdups::template_ends pair with both ends unmapped", "[template_ends]"){
  // Confirm that ends of pair with both ends unmapped are calculated as expected
  SamRecord pair1_r1("HWI-ST1213:151:C1DTBACXX:2:1112:10072:37225\t77\t*\t0\t0\t*\t*\t0\t0\tAGCTACCGCCTGCCCCAGGACCTGACCCTCTTGAAATCCAACCACTAACGTCCCCCACCACGTCTCAACTTAGCAGTTCGGCACCAAGCTGCGCACAAACT\tCCCFFFFFHHHHHJJJJJJJIJJJJJJJIJJJJJIJJJJJJJJJJJJJJJGIJJJHHFFFFCADEDEDDDDCCDDCCCCDDDDDDDADDDCDDDDDDDBD4\tPG:Z:MarkDuplicates\tAS:i:0\tXS:i:0");
  SamRecord pair1_r2("HWI-ST1213:151:C1DTBACXX:2:1112:10072:37225\t141\t*\t0\t0\t*\t*\t0\t0\tCTGAGCCCAAGTGGGCTTATATGGTGTTATTTCCACTTGAGGTGTCTGATGGGGAGTGGGTCCCACTGGATGATTAAATGGGATTGGCTGCATCCGCAAGG\tCCCFFFFFHHHDFIJJJJIJJJJJCHHHJJJJJJJJJJJEHIBFFHIJIJHHJJGICHIIEHHIGHHHHFFFFFFEEEEEEDCCDDDDDDBCDCCDDBBBB\tPG:Z:MarkDuplicates\tAS:i:0\tXS:i:0");
  std::vector<SamRecord> qn1 { pair1_r1, pair1_r2 };
  auto ends1 { ends_str(template_ends(qn1), 2) };
  CHECK(ends1 == unmapped + "_" + unmapped);
}

TEST_CASE("markdups::template_ends single FF", "[template_ends]"){
  // Confirm that calculated ends of aligned single-end reads with same RNAME,
  // POS and orientation, are the same.
  SamRecord sr1("NB551151:333:fake1:4:21501:17121:9587\t0\tchr1\t4498454\t255\t73M\t*\t0\t0\tGTCTTCCGTTCACTACACCTTTCAATCCTGGATCACAGGGCTTTCCAGCCTTGACTACATACTTACGAATAAT\tAAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEE\tNH:i:1\tHI:i:1\tAS:i:71\tnM:i:0\tNM:i:0\tMD:Z:73\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e");
  SamRecord sr2("NB551151:333:fake1:4:21506:21331:18924\t0\tchr1\t4498454\t255\t73M\t*\t0\t0\tGTCTTCCGTTCACTACACCTTTCAATCCTGGATCACAGGGCTTTCCAGCCTTGACTACATACTTACGAATAAT\t/EAEEEAEAEEEEEEEEEEEEEEEEAEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAA\tNH:i:1\tHI:i:1\tAS:i:74\tnM:i:0\tNM:i:0\tMD:Z:76\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e");
  std::vector<SamRecord> qn1 { sr1 }, qn2 { sr2 };
  auto ends1 { ends_str(template_ends(qn1), 1) },
       ends2 { ends_str(template_ends(qn2), 1) };
  CHECK(ends1 == "chr1F4498454");
  CHECK(ends1 == ends2);
}

TEST_CASE("markdups::template_ends single FR", "[template_ends]"){
  // Confirm that calculated ends of aligned single-end reads with same RNAME,
  // POS and opposite orientation, are the different.
  SamRecord sr1("NB551151:333:fake1:4:21501:17121:9587\t0\tchr1\t4498454\t255\t73M\t*\t0\t0\tGTCTTCCGTTCACTACACCTTTCAATCCTGGATCACAGGGCTTTCCAGCCTTGACTACATACTTACGAATAAT\tAAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEE\tNH:i:1\tHI:i:1\tAS:i:71\tnM:i:0\tNM:i:0\tMD:Z:73\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e");
  SamRecord sr2("NB551151:333:fake1:4:21506:21331:18924\t16\tchr1\t4498454\t255\t73M\t*\t0\t0\tGTCTTCCGTTCACTACACCTTTCAATCCTGGATCACAGGGCTTTCCAGCCTTGACTACATACTTACGAATAAT\t/EAEEEAEAEEEEEEEEEEEEEEEEAEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAA\tNH:i:1\tHI:i:1\tAS:i:74\tnM:i:0\tNM:i:0\tMD:Z:76\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e");
  std::vector<SamRecord> qn1 { sr1 }, qn2 { sr2 };
  auto ends1 { ends_str(template_ends(qn1), 1) },
       ends2 { ends_str(template_ends(qn2), 1) };
  CHECK(ends1 == "chr1F4498454");
  CHECK(ends1 != ends2);
}

TEST_CASE("markdups::template_ends paired FR FR", "[template_ends]"){
  // Confirm that calculated ends of a duplicate pair in same orientation are
  // the same.
  SamRecord pair1_r1("HWI-ST1213:151:C1DTBACXX:2:1309:19267:90133\t99\tchr1\t13583\t22\t101M\t=\t13759\t277\tGAAACCCAGGAGAGTGTGGAGTCCAGAGTGATGCCAGGACCCAGGCACAGGCATTAGTGCCCGTTGGAGAAAACAGGGGAATCCCGAAGAAATGGTGGGTC\tCCCFFFFFHHHHHJEGGIIJJHIHIJJJGHIJJJIJJJJJJJJIJJIJIJIJJJGJIGGIJJHHGFFDDDDEEEDDDDDBBDDCDDDDDDCDCDD@CDD53\tMD:Z:30T70\tPG:Z:MarkDuplicates\tDI:i:1\tNM:i:1\tAS:i:96\tDS:i:2\tXS:i:96");
  SamRecord pair1_r2("HWI-ST1213:151:C1DTBACXX:2:1309:19267:90133\t147\tchr1\t13759\t22\t101M\t=\t13583\t-277\tGGGCCAGGCTTCTCACTGGGCCTCTGCAGAAGGCTGCCATTTGTCCTGCCCACCGTCTTAGAAGCGAGACGGAGCAGACTCATCTGCTACTGCCCTTTCTA\t?B?DBDBCCDDCDDCDBDDDDDDDDDDDDDEDDDFFFFFFHHHHJIIGGHFAJJJJJJJJIJIIJJJIHJJJIJJJJIJIJJJJJJJIHHGGHFFFFFCCC\tMD:Z:29G24T24C21\tPG:Z:MarkDuplicates\tDI:i:1\tNM:i:3\tAS:i:86\tDS:i:2\tXS:i:86");
  SamRecord pair2_r1("HWI-ST1213:151:C1DTBACXX:2:1108:11333:45933\t1123\tchr1\t13583\t22\t101M\t=\t13759\t277\tGAAACCCAGGAGAGTGTGGAGTCCAGAGTGATGCCAGGACCCAGGCACAGGCATTAGTGCCCGTTGGAGAAAACAGGGGAATCCCGAAGAAATGGTGGGTC\tCCCFFFFFHHHGHJFHFHIIJGIIJIJJFGIJJIJJJJIIJJIJIIJIIJJJEGIJJIIDHGFHEFEEFDEEEEDD??BDDDDDDDDDDBDDDDCCDDD##\tMD:Z:30T70\tPG:Z:MarkDuplicates\tDI:i:1\tNM:i:1\tAS:i:96\tDS:i:2\tXS:i:96");
  SamRecord pair2_r2("HWI-ST1213:151:C1DTBACXX:2:1108:11333:45933\t1171\tchr1\t13759\t22\t101M\t=\t13583\t-277\tGGGCCAGGCTTCTCACTGGGCCTCTGCAGAAGGCTGCCATTTGTCCTGCCCACCGTCTTAGAAGCGAGACGGAGCAGACTCATCTGCTACTGCCCTTTCTA\tDDDDBCBC?C>CACCB?DDBDBCC>CCDDDEC@CFFFFDFGHHHIIIGGF=7HEHGIIGGGIIIJJJJIIIJJIGHHEGJJJJGJIIJHHHHGFFFFFCC@\tMD:Z:29G24T24C21\tPG:Z:MarkDuplicates\tDI:i:1\tNM:i:3\tAS:i:86\tDS:i:2\tXS:i:86");
  std::vector<SamRecord> qn1 { pair1_r1, pair1_r2 }, qn2 { pair2_r1, pair2_r2 };
  auto ends1 { ends_str(template_ends(qn1), 2) },
       ends2 { ends_str(template_ends(qn2), 2) };
  CHECK(ends1 == "chr1F13583_chr1R13860");
  CHECK(ends1 == ends2);
}

TEST_CASE("markdups::template_ends paired FR RF", "[template_ends]"){
  // Confirm that calculated ends of a duplicate pair in opposite orientation
  // are the same.
  SamRecord pair1_r1("HWI-ST1213:151:C1DTBACXX:2:1207:6478:58471\t99\tchr1\t564691\t0\t101M\t=\t564887\t297\tCGATTCCGCTACGACCAACTCATACACCTCCTATGAAAAAACTTCCTACCACTCACCCTAGCATTACTTATATGATATGTCTCCATACCCATTACAATCTC\tCCCFFFFFHHHHHJJIJIJIIJJIJJJJJJJJIIIJIIJGIIJIJJIJJFIJIJJHHHEFFFFFFCEDCEDEEFDCCFDDDDDDDACDDDDCDDDD>CCDC\tXA:Z:chrMT,+4141,101M,0;\tMD:Z:101\tPG:Z:MarkDuplicates\tDI:i:10\tNM:i:0\tAS:i:101\tDS:i:2\tXS:i:101");
  SamRecord pair1_r2("HWI-ST1213:151:C1DTBACXX:2:1207:6478:58471\t147\tchr1\t564887\t0\t101M\t=\t564691\t-297\tATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACC\tCDCA@@C<DB?>D@DC@BDCEDDDDDDBDDDDDDDDADFFFHHE?HIGHE@@GBJJIIHGJJIJJJJIJJJIJIIJIIJJJIJJJIJIHGHHHFFFDD=CB\tXA:Z:chrMT,-4337,101M,0;\tMD:Z:101\tPG:Z:MarkDuplicates\tDI:i:10\tNM:i:0\tAS:i:101\tDS:i:2\tXS:i:101");
  SamRecord pair2_r1("HWI-ST1213:151:C1DTBACXX:2:2105:17263:64957\t1107\tchr1\t564887\t0\t101M\t=\t564691\t-297\tATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACC\t>DCA:4(85229?<83<@C@;;C@55:;3AA>5;,9/C?@C>C?7?:C<H@5A7FDCCGACHAB9D9CGEGDBDFD<@EIHHCGDC??:AAFFDDBDD?@@\tXA:Z:chrMT,-4337,101M,0;\tMD:Z:101\tPG:Z:MarkDuplicates\tDI:i:10\tNM:i:0\tAS:i:101\tDS:i:2\tXS:i:101");
  SamRecord pair2_r2("HWI-ST1213:151:C1DTBACXX:2:2105:17263:64957\t1187\tchr1\t564691\t0\t101M\t=\t564887\t297\tCGATTCCGCTACGACCAACTCATACACCTCCTATGAAAAAACTTCCTACCACTCACCCTAGCATTACTTATATGATATGTCTCCATACCCATTACAATCTC\t?1=BD;;3C::6DF:1EGGBHII<<FHGDGBFH4DGDAFHGIGIADA9@@CD3C3?>EE@9BDD>A>CCCEACCC@D;5>;>>CCDCC?CCDDCCC>@>AC\tXA:Z:chrMT,+4141,101M,0;\tMD:Z:101\tPG:Z:MarkDuplicates\tDI:i:10\tNM:i:0\tAS:i:101\tDS:i:2\tXS:i:101");
  std::vector<SamRecord> qn1 { pair1_r1, pair1_r2 }, qn2 { pair2_r1, pair2_r2 };
  auto ends1 { ends_str(template_ends(qn1), 2) },
       ends2 { ends_str(template_ends(qn2), 2) };
  CHECK(ends1 == "chr1F564691_chr1R564988");
  CHECK(ends1 == ends2);
}

TEST_CASE("markdups::template_ends pairs with soft-clipping", "[template_ends]"){
  // Confirm that soft clipping at fwd and reverse ends of reads in a duplicate
  // pair is accounted for
  SamRecord pair1_r1("HWI-ST1213:151:C1DTBACXX:2:1202:18193:22243\t83\tchr1\t725905\t19\t101M\t=\t725717\t-289\tACGAATGGAATGGAATGGAATGGACTCGAATGGAATTGAATGCAATGGAATGGACCTGAGAGGATTGGAATGGAATTGACTTGAATGGAATTGAATGGAAC\t=HHEHGGJHEIJJJJIIJJIJIHHIHIGHJJJJJJJJJJJIJJJIJJJJJJIHGIJIJJJJJJJJJJJJJJJJJJJIJJJJJJJJJJJHHHHHFFFFFCCC\tMD:Z:0T35G2C15T0C0A1A0T3A10C5C18T0\tPG:Z:MarkDuplicates\tDI:i:44\tNM:i:12\tAS:i:49\tDS:i:2\tXS:i:44");
  SamRecord pair1_r2("HWI-ST1213:151:C1DTBACXX:2:1202:18193:22243\t163\tchr1\t725717\t19\t6S39M10I46M\t=\t725905\t289\tGAACTCGAATGTAATACAATGGAATTTAATGGAATGGATTCTAATGGAATAGAAAGGAATGGACTCGAATGGAATAGAATGGAATGGACTCGAATGGAATG\tC@CFFFFFHHHFHJJIJJJJJJIJJJJJIJJJJJJJJIIJJJJIJJJGIEAFGIIJJJJJJJJJJJJIJJJJJJJHGHHDHC>DEFFDEDDDDDDDDCCD@\tMD:Z:5G26C16T18G0A15\tPG:Z:MarkDuplicates\tDI:i:44\tNM:i:15\tAS:i:44\tDS:i:2\tXS:i:55");
  SamRecord pair2_r1("HWI-ST1213:151:C1DTBACXX:2:1107:3646:4446\t1123\tchr1\t725721\t30\t10S32M10I34M15S\t=\t725895\t275\tGAACACGAAAGGAATACAAGGGAATTTAATGGAATGGACTCTAATGGAATGAAATGGAATGGACTTGAATGGAATATAATGGAAGATATTAGAATGGAATA\tCCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJIJJJJIIJJJIJIJJIIJIJJJJJJEHHIJIJJJJIGHHHHHFFDEFFFEEEECDDDFEEDDDDDCDD>\tMD:Z:9T46G9\tPG:Z:MarkDuplicates\tDI:i:44\tNM:i:12\tAS:i:40\tDS:i:2\tXS:i:43");
  SamRecord pair2_r2("HWI-ST1213:151:C1DTBACXX:2:1107:3646:4446\t1171\tchr1\t725895\t30\t40M10D51M10S\t=\t725721\t-275\tTGGAATGGACTCGAATGGAATGGTATGGAATGGACTCGAATGCAATGGAATGTACTCAAATGGAATGCTATGGAATTGACTCGAGTGGAATGGAATAGAAT\t;FFFFEEHFHHIIJJIJIJJIIIJIJJJJJIJJIJJIJJJJJJJJJIGIGGIGDJJGJJJIJJJJJJJJJJJJJJJJJHJJJJJJJJJHHHHHFFFFFCCC\tMD:Z:8T14A16^TGGAATGGAC12G14G0A6C8A6\tPG:Z:MarkDuplicates\tDI:i:44\tNM:i:17\tAS:i:40\tDS:i:2\tXS:i:35");
  std::vector<SamRecord> qn1 { pair1_r1, pair1_r2 }, qn2 { pair2_r1, pair2_r2 };
  auto ends1 { ends_str(template_ends(qn1), 2) },
       ends2 { ends_str(template_ends(qn2), 2) };
  CHECK(ends1 == "chr1F725711_chr1R726006");
  CHECK(ends1 == ends2);
}

TEST_CASE("markdups::process_qname_group not single", "[process_qname_group]"){
  // Confirm that exception is thrown if reads_per_template is 1 but qname
  // group does not have exactly 1 primary alignment
  SamRecord pair1_r1("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25");
  SamRecord pair1_r2("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t147\tchr1\t93578228\t60\t101M\t=\t93578030\t-299\tAACAACAACAAAAAATTTGGTATTTCTAAGATGAAATGGCCAAGGCTTTCTAGTCAATTGGATTTAGAGTAAAGGAGACTATAGAAGATTACTAAGCTATA\tBDDDDDDEDDFHHHHHHJIJIJJJJJIJJJIJJJJJJJJJIIJJIJJJIIIJJIJJJJJIJJJJIJJJJJJJJJJJJJJJIJJJJJJJHHHHHFFFFFCCB\tNM:i:0\tAS:i:101\tXS:i:21");
  std::vector<SamRecord> qn1 { pair1_r1, pair1_r2 };
  std::stringstream os;
  bloomfilter::BloomFilter bf(0.001, 1000);
  uint64_t n_tpl_dup {0}, n_aln_dup {0};
  size_t reads_per_template = 1;
  CHECK_THROWS_WITH(
      process_qname_group(qn1, os, bf, n_tpl_dup, n_aln_dup, reads_per_template),
      Catch::Matchers::Contains("Input is not single reads?"));
}

TEST_CASE("markdups::process_qname_group not paired", "[process_qname_group]"){
  // Confirm that exception is thrown if reads_per_template is 2 but qname
  // group does not have exactly 2 primary alignments
  SamRecord sr1("HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tAS:i:101\tXS:i:25");
  std::vector<SamRecord> qn1 { sr1 };
  std::stringstream os;
  bloomfilter::BloomFilter bf(0.001, 1000);
  uint64_t n_tpl_dup {0}, n_aln_dup {0};
  size_t reads_per_template = 2;
  CHECK_THROWS_WITH(
      process_qname_group(qn1, os, bf, n_tpl_dup, n_aln_dup, reads_per_template),
      Catch::Matchers::Contains("Input is not paired or not qname-grouped?"));
}

TEST_CASE("markdups::process_qname_group strip_previous==false", "[process_qname_group]"){
  // Confirm that when strip_previous is false (default) that previously marked
  // duplicate flag remains on read no longer considered duplicate
  SamRecord sr1("NB551151:333:fake1:4:21501:17121:9587\t1024\tchr1\t4498454\t255\t73M\t*\t0\t0\tGTCTTCCGTTCACTACACCTTTCAATCCTGGATCACAGGGCTTTCCAGCCTTGACTACATACTTACGAATAAT\tAAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEE\tNH:i:1\tHI:i:1\tAS:i:71\tnM:i:0\tNM:i:0\tMD:Z:73\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e\tPG:Z:MarkDuplicates");
  SamRecord sr2("NB551151:333:fake1:4:21506:21331:18924\t0\tchr1\t4498454\t255\t73M\t*\t0\t0\tGTCTTCCGTTCACTACACCTTTCAATCCTGGATCACAGGGCTTTCCAGCCTTGACTACATACTTACGAATAAT\t/EAEEEAEAEEEEEEEEEEEEEEEEAEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAA\tNH:i:1\tHI:i:1\tAS:i:74\tnM:i:0\tNM:i:0\tMD:Z:76\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e");
  std::vector<SamRecord> qn1 { sr1 };
  std::vector<SamRecord> qn2 { sr2 };
  REQUIRE(sr1.flag() == 1024);
  REQUIRE_THAT(sr1.buffer, Catch::Matchers::Contains("PG:Z:MarkDuplicates"));
  REQUIRE(sr2.flag() == 0); 
  REQUIRE_THAT(sr2.buffer, !Catch::Matchers::Contains("PG:Z:"));
  std::stringstream os;
  bloomfilter::BloomFilter bf(0.001, 1000);
  uint64_t n_tpl_dup {0}, n_aln_dup {0};
  process_qname_group(qn1, os, bf, n_tpl_dup, n_aln_dup, 1, false);
  process_qname_group(qn2, os, bf, n_tpl_dup, n_aln_dup, 1, false);
  qn1[0].parse();
  qn2[0].parse();
  // original dupe flag and PG remain on first read even though it was processed first
  CHECK(qn1[0].flag() == 1024);
  CHECK_THAT(qn1[0].buffer, Catch::Matchers::Contains("PG:Z:MarkDuplicates"));
  CHECK(qn2[0].flag() == 1024);
  CHECK_THAT(qn2[0].buffer, Catch::Matchers::Contains("PG:Z:streammd"));
}

TEST_CASE("markdups::process_qname_group strip_previous==true", "[process_qname_group]"){
  // Confirm that when strip_previous is true that previously marked duplicate
  // flag is removed from read no longer considered duplicate.
  SamRecord sr1("NB551151:333:fake1:4:21501:17121:9587\t1024\tchr1\t4498454\t255\t73M\t*\t0\t0\tGTCTTCCGTTCACTACACCTTTCAATCCTGGATCACAGGGCTTTCCAGCCTTGACTACATACTTACGAATAAT\tAAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEE\tNH:i:1\tHI:i:1\tAS:i:71\tnM:i:0\tNM:i:0\tMD:Z:73\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e\tPG:Z:MarkDuplicates");
  SamRecord sr2("NB551151:333:fake1:4:21506:21331:18924\t0\tchr1\t4498454\t255\t73M\t*\t0\t0\tGTCTTCCGTTCACTACACCTTTCAATCCTGGATCACAGGGCTTTCCAGCCTTGACTACATACTTACGAATAAT\t/EAEEEAEAEEEEEEEEEEEEEEEEAEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAA\tNH:i:1\tHI:i:1\tAS:i:74\tnM:i:0\tNM:i:0\tMD:Z:76\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e");
  std::vector<SamRecord> qn1 { sr1 };
  std::vector<SamRecord> qn2 { sr2 };
  REQUIRE(sr1.flag() == 1024);
  REQUIRE_THAT(sr1.buffer, Catch::Matchers::Contains("PG:Z:MarkDuplicates"));
  REQUIRE(sr2.flag() == 0); 
  REQUIRE_THAT(sr2.buffer, !Catch::Matchers::Contains("PG:Z:"));
  std::stringstream os;
  bloomfilter::BloomFilter bf(0.001, 1000);
  uint64_t n_tpl_dup {0}, n_aln_dup {0};
  process_qname_group(qn1, os, bf, n_tpl_dup, n_aln_dup, 1, true);
  process_qname_group(qn2, os, bf, n_tpl_dup, n_aln_dup, 1, true);
  qn1[0].parse();
  qn2[0].parse();
  // original dupe flag and PG stripped from first read
  CHECK(qn1[0].flag() == 0);
  CHECK_THAT(qn1[0].buffer, !Catch::Matchers::Contains("PG:Z:MarkDuplicates"));
  CHECK(qn2[0].flag() == 1024);
  CHECK_THAT(qn2[0].buffer, Catch::Matchers::Contains("PG:Z:streammd"));
}

TEST_CASE("markdups::process_input_stream single reads", "[process_input_stream]"){
  // Confirm that duplicates are marked as expected on single-end reads.
  std::stringstream in {
    "NB551151:333:fake1:4:21501:17121:9587\t0\tchr1\t4498454\t255\t73M\t*\t0\t0\tGTCTTCCGTTCACTACACCTTTCAATCCTGGATCACAGGGCTTTCCAGCCTTGACTACATACTTACGAATAAT\tAAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAEE\tNH:i:1\tHI:i:1\tAS:i:71\tnM:i:0\tNM:i:0\tMD:Z:73\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e\n"
    "NB551151:333:fake1:4:21506:21331:18924\t0\tchr1\t4498454\t255\t73M\t*\t0\t0\tGTCTTCCGTTCACTACACCTTTCAATCCTGGATCACAGGGCTTTCCAGCCTTGACTACATACTTACGAATAAT\t/EAEEEAEAEEEEEEEEEEEEEEEEAEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEAA\tNH:i:1\tHI:i:1\tAS:i:74\tnM:i:0\tNM:i:0\tMD:Z:76\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e\n"
  "NB551151:333:fake1:3:13509:15168:5615\t16\tchr1\t4502637\t255\t75M1S\t*\t0\t0\tGGATGTCCTCACCCTAGTCAGACCAATGTTGGCCCACACATTCTTTGCAGGAACCCCACAGAAATTGTAGTACCCC\tEEEEEEEEEAEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEAEEEEEEEEEAAAAA\tNH:i:1\tHI:i:1\tAS:i:73\tnM:i:0\tNM:i:0\tMD:Z:75\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ba5a08e7-a5f7-48ce-b8f8-c0f223463c29\n"
  "NB551151:333:fake1:4:12509:22440:18769\t16\tchr1\t4505378\t255\t76M\t*\t0\t0\tAAATCAACCCTTTCCTCTTCAACTTGTCTTGGTCATGGTGTTTCATCACAGCAATAGTAACCCTAACTAAAACAGG\tEEE<E/AEEAAEE/EEEEEA6/EEEEEEEE6AAEEEAEEEEE/EEEAEEEEEEEEEEEE//EAEEAEEEA/AA/AA\tNH:i:1\tHI:i:1\tAS:i:74\tnM:i:0\tNM:i:0\tMD:Z:76\tjM:B:c,-1\tjI:B:i,-1\tRG:Z:ec6fd08a-5874-476c-a6c6-24ba72aa308e\n" };
  std::stringstream out;
  bloomfilter::BloomFilter bf(0.001, 1000);
  std::vector<std::string> cli_args { "dummy", "args" };
  auto result = process_input_stream(in, out, bf, cli_args, 1);
  CHECK(result.templates == 4);
  CHECK(result.templates_marked_duplicate == 1); // the second read
  CHECK(result.alignments == 4);
  CHECK(result.alignments_marked_duplicate == 1);
}

TEST_CASE("markdups::process_input_stream paired reads", "[process_input_stream]"){
  // Confirm that duplicates are marked as expected on paired-end reads.
  std::stringstream in {
    "HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t99\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHHJJJJJGIJJJJJIHHIJJJJJJIJJJJJJJJJIJJJJJHIIJGHFHHGIIJJGJIJIJJJEHHEHEHDDEFFFEEEDEEEDDDEEDCC\tNM:i:0\tMD:Z:101\tAS:i:101\tXS:i:25\n"
    "HWI-ST1213:151:C1DTBACXX:2:1101:2189:99680\t147\tchr1\t93578228\t60\t101M\t=\t93578030\t-299\tAACAACAACAAAAAATTTGGTATTTCTAAGATGAAATGGCCAAGGCTTTCTAGTCAATTGGATTTAGAGTAAAGGAGACTATAGAAGATTACTAAGCTATA\tBDDDDDDEDDFHHHHHHJIJIJJJJJIJJJIJJJJJJJJJIIJJIJJJIIIJJIJJJJJIJJJJIJJJJJJJJJJJJJJJIJJJJJJJHHHHHFFFFFCCB\tNM:i:0\tMD:Z:101\tAS:i:101\tXS:i:21\n"
    "HWI-ST1213:151:C1DTBACXX:2:2106:1715:59565\t83\tchr1\t93578228\t60\t101M\t=\t93578030\t-299\tAACAACAACAAAAAATTTGGTATTTCTAAGATGAAATGGCCAAGGCTTTCTAGTCAATTGGATTTAGAGTAAAGGAGACTATAGAAGATTACTAAGCTATA\t?DDDDDDEDDFGHGFHHGIJJJIJJIIJJIIJIJJJIJJIIIJJJIIJIJJJJIJIJJIJJJIIJJJJJJJJJJJJJIJJJJJJJJJJHHHGHFFFFFCCC\tNM:i:0\tMD:Z:101\tAS:i:101\tXS:i:21\n"
    "HWI-ST1213:151:C1DTBACXX:2:2106:1715:59565\t163\tchr1\t93578030\t60\t101M\t=\t93578228\t299\tCATCTAATGCTGTTTTGGTTTCTGTGAAATGATATACTCTTGCATTGCTGGTGGCAGTGTAAATTTCTATTTTGGTGGTTTAGCATTTGCATTAAATGCAG\tCCCFFFFFHHHHFHIJJJFHIJJJJJJJJJJJJJJJJJJJJIJJJIIJIJJHHIIJJHIGIIGIJIJJIJJJIIJHIHHFFFFFFFEEEDEDEDDDDDDDC\tNM:i:0\tMD:Z:101\tAS:i:101\tXS:i:25\n"
  };
  std::stringstream out;
  bloomfilter::BloomFilter bf(0.001, 1000);
  std::vector<std::string> cli_args { "dummy", "args" };
  auto result = process_input_stream(in, out, bf, cli_args, 2);
  CHECK(result.templates == 2);
  CHECK(result.templates_marked_duplicate == 1); // the second pair 
  CHECK(result.alignments == 4);
  CHECK(result.alignments_marked_duplicate == 2);
}

TEST_CASE("markdups::process_input_stream unmapped", "[process_input_stream]"){
  // Confirm that paired records with both read and mate unmapped still appear
  // in the output
  std::stringstream in {
    "HWI-ST1213:151:C1DTBACXX:2:2207:13476:31678\t77\t*\t0\t0\t*\t*\t0\t0\tGGCCGACAACCAGATTATGGGGCCCTCGAGCTATGCTTCTTTTGTGGTACGGGGGGAGAACCTGGTCACTGCCGTGAGCTACGGGCGCGTGATGCGTACGT\tCCCFFFFFHHHHHJJJIJIJJJIJJJJJJJJJIJJJJJJJJIJJGHJB@DGEHFDDDDDBDDDDDDDCDDDDDDB@DDCDDDDBDB@B>9B?CDCD88>BC\tAS:i:0\tXS:i:0\n"
    "HWI-ST1213:151:C1DTBACXX:2:2207:13476:31678\t141\t*\t0\t0\t*\t*\t0\t0\tCTTCTGTTGAGGGGGTATGGGGACTGAGTGTCATTGTACATCTTTTGCAGGCTTTCCACGGCCACCGCGTGGTTGCCCAGCTTGATGACGGCGGCTG\tCCCFFFFFGHHHHJJ?FHIIJJHIJJIJHIIIJJJJGIJIJJJJJJIJIJGHIJIJHHHHFFDDDDDDDBDDADDDDDDDDDDDDDDEDDBDD><.9\tAS:i:0\tXS:i:0\n"
  };
  std::stringstream out;
  bloomfilter::BloomFilter bf(0.001, 1000);
  std::vector<std::string> cli_args { "dummy", "args" };
  auto result = process_input_stream(in, out, bf, cli_args, 2);
  CHECK(result.templates == 1);
  CHECK(result.templates_marked_duplicate == 0); 
  CHECK(result.alignments == 2);
  CHECK(result.alignments_marked_duplicate == 0);
}

/*
   Confirm that process_input_stream operates as expected on a larger input
   file: 4058 alignments from 2027 templates, with a variety of different
   flags, orientations, and complex CIGAR strings.
 
   The test reference output was generated from the test input by Picard
   MarkDuplicates (2.23.8).
 
   The input SAM records in 'test.paired_full.sam' have been ordered with
   higher quality reads in a duplicate set occuring first so that in this
   case we expect the output of streammd to EXACTLY match that from
   Picard MarkDuplicates (assuming no false positives from streammd,
   which happens to be the case for this data with default --fp-rate and
   seeds).
 
   In the general case we expect only that the duplicate counts should
   match, since Picard MarkDuplicates picks the highest quality read in a
   set as the original, where streammd must pick the first.
 */
TEST_CASE("markdups::process_input_stream full SAM", "[process_input_stream]") {
  std::ifstream testinstrm, expectinstrm;
  testinstrm.open("resources/test.paired_full.sam");
  expectinstrm.open("resources/test.paired_full.picardmd.sam");
  std::ostringstream testoutstrm;
  std::map<
    std::tuple<std::string, std::string, size_t>,
    uint16_t> expected_flags, marked_flags;
  for (SamRecord sr; std::getline(expectinstrm, sr.buffer); ) {
    if (sr.buffer[0] != '@') {
      sr.parse();
      expected_flags[
        std::make_tuple(sr.qname(), sr.rname(), sr.pos())] = sr.flag();
    }
  }
  bloomfilter::BloomFilter bf(0.000001, 1000000);
  std::vector<std::string> cli_args { "dummy", "args" };
  auto result = process_input_stream(testinstrm, testoutstrm, bf, cli_args, 2);
  CHECK(result.templates == 2027);
  CHECK(result.templates_marked_duplicate == 1018); 
  CHECK(result.alignments == 4058);
  CHECK(result.alignments_marked_duplicate == 2039);
  CHECK_THAT(float(result.templates_marked_duplicate)/
             result.templates,
             Catch::Matchers::WithinAbs(0.5022, 0.0001));
  auto outlines { std::istringstream(testoutstrm.str()) };
  for (SamRecord sr; std::getline(outlines, sr.buffer); ) {
    if (sr.buffer[0] != '@') {
      sr.parse();
      marked_flags[
        std::make_tuple(sr.qname(), sr.rname(), sr.pos())] = sr.flag();
    }
  }
  // check that all alignments have expected flag
  for (auto const& [key, val]: expected_flags) {
    CHECK(marked_flags[key] == val);
  }
}
