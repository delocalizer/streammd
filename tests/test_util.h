#include <iostream>
#include <map>
#include "markdups.h"

#include <catch2/catch.hpp>

// (QNAME, RNAME, POS)
typedef std::tuple<std::string, std::string, size_t> SamRecordId;
typedef std::map<SamRecordId, uint16_t> SamRecordFlags;

using namespace markdups;

// Populate SAM flags map from SAM record input stream
void record_flags(std::istream& instream, SamRecordFlags& flags);

// Run process_input_stream against input SAM file; check SAM output against the
// reference file. The metrics struct is returned for possible further checks.
metrics test_streammd(
    std::string input_path,
    std::string reference_path,
    size_t reads_per_template=2,
    double p=0.000001,
    uint64_t=1000000
    );
