#include <fstream>
#include <iostream>

#include <catch2/catch.hpp>

#include "markdups.h"
#include "test_util.h"


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

// Run process_input_stream against input SAM file; check SAM output against the
// reference file. The metrics struct is returned for possible further checks.
metrics test_streammd(
    std::string input_path,
    std::string reference_path,
    size_t reads_per_template,
    double p,
    uint64_t n
    ){
  std::ifstream test_input, expected_output;
  SamRecordFlags marked_flags, expected_flags;
  test_input.open(input_path);
  expected_output.open(reference_path);
  record_flags(expected_output, expected_flags);

  std::ostringstream test_output;
  bloomfilter::BloomFilter bf(p, n);
  std::vector<std::string> cli_args { "dummy", "args" };
  auto result = process_input_stream(
      test_input, test_output, bf, cli_args, reads_per_template);
  auto outlines { std::istringstream(test_output.str()) };
  record_flags(outlines, marked_flags);

  for (auto const& [key, val]: expected_flags) {
    INFO(fmt::format("SamRecordId: {}:{}:{}",
                     std::get<0>(key), std::get<1>(key), std::get<2>(key)));
    CHECK(marked_flags[key] == val);
  }
  return result;
}
