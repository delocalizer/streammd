#ifndef STREAMMD_MARKDUPS_H_
#define STREAMMD_MARKDUPS_H_

#include <cstdint>
#include <regex>
#include <string>
#include <tuple>

#include "bloomfilter.h"

namespace markdups {

// DEL sorts last in ASCII
const char DEL { 127 };
const char SAM_delimiter { '\t' };
const float default_p { 0.000001 };
const std::regex re_cigar { R"((?:(\d+)([MIDNSHPX=])))" };
const std::regex re_leading_s { R"(^(\d+)S)" };                       
const std::regex re_trailing_s { R"((\d+)S$)" };
const std::set consumes_reference { 'M', 'D', 'N', '=', 'X' };
const std::string default_metrics { "streammd-metrics.json" };
const std::tuple<std::string, uint32_t, char> unmapped { std::string(1, DEL), 0, DEL };
const uint32_t log_interval { 1000000 };
const uint64_t default_n { 1000000000 };
const unsigned short flag_unmapped = 4;
const unsigned short flag_reverse = 16;
const unsigned short flag_secondary = 256;
const unsigned short flag_duplicate = 1024;
const unsigned short flag_supplementary = 2048;

inline std::string join(
    const std::vector<std::string>& elems,
    const char& sep,
    int reserve=1024) {
  std::string joined;
  joined.reserve(reserve);
  unsigned long i {1}, imax { elems.size() };
  for (auto & elem : elems) {
    joined += elem;
    if (i < imax) { joined += sep; }
    i++;
  }
  return joined;
}

inline std::vector<std::string> split(
    const std::string& joined, const char& sep) {
  std::vector<std::string> elems;
  std::stringstream joined_{ joined };
  std::string tkn;
  while(std::getline(joined_, tkn, sep)) {
    elems.push_back(tkn);
  }
  return elems;
}

void mark_duplicates(
    std::vector<std::vector<std::string>>& qname_group,
    int& reads_per_template,
    bloomfilter::BloomFilter& bf);

void readends(
    std::vector<std::vector<std::string>>& qname_group,
    std::vector<std::tuple<std::string, uint32_t, char>>& ends);

void process(
    std::istream& in,
    std::ostream& out,
    int reads_per_template,
    bloomfilter::BloomFilter& bf);

void write(
    std::vector<std::vector<std::string>>& qname_group,
    std::ostream& out);
}
#endif // STREAMMD_MARKDUPS_H_
  
