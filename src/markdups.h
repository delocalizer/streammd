#ifndef STREAMMD_MARKDUPS_H_
#define STREAMMD_MARKDUPS_H_

#include <cstdint>
#include <regex>
#include <string>
#include <tuple>

#include "bloomfilter.h"

using end_t = std::tuple<std::string, uint32_t, char>;

namespace markdups {

// DEL sorts last in ASCII
const char DEL { 127 };
const char SAM_delimiter { '\t' };
const float default_p { 0.000001 };
const std::regex re_cigar { R"((?:(\d+)([MIDNSHPX=])))" };
const std::regex re_leading_s { R"(^(\d+)S)" };                       
const std::regex re_trailing_s { R"((\d+)S$)" };
const std::set<std::string> consumes_reference { "M", "D", "N", "=", "X" };
const std::string pgid { "streammd"  };
const std::string pgtag { "PG:Z" };
const std::string pgtag_val { pgtag + ":" + pgid };
const std::string default_metrics { pgid + "-metrics.json" };
const end_t unmapped { std::string(1, DEL), -1, DEL };
const uint32_t log_interval { 1000000 };
const uint64_t default_n { 1000000000 };
const unsigned short flag_unmapped = 4;
const unsigned short flag_reverse = 16;
const unsigned short flag_secondary = 256;
const unsigned short flag_duplicate = 1024;
const unsigned short flag_supplementary = 2048;
const unsigned short sam_opts_idx = 11;

inline std::string join(
    const std::vector<std::string>& elems,
    const char sep,
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
    const std::string& joined, const char sep) {
  std::vector<std::string> elems;
  std::stringstream joined_{ joined };
  std::string tkn;
  while(std::getline(joined_, tkn, sep)) {
    elems.push_back(tkn);
  }
  return elems;
}

inline std::string ends_to_string(std::vector<end_t> ends) {
  std::string ends_str;
  for (auto end : ends) {
    ends_str += std::get<0>(end);
    // Place F|R next because rname can end in a digit and we need string
    // values that distinguish between ("chr1", 1234) and ("chr11", 234) 
    ends_str += std::get<2>(end);
    ends_str += std::to_string(std::get<1>(end));
    // Ditto, in case the next rname starts with a digit
    ends_str += "_";
  }
  return ends_str;
}

void process_input_stream(
    std::istream& in,
    std::ostream& out,
    bloomfilter::BloomFilter& bf,
    unsigned reads_per_template = 2,
    bool strip_previous = false);

void process_qname_group(
    std::vector<std::vector<std::string>>& qname_group,
    std::ostream& out,
    bloomfilter::BloomFilter& bf,
    unsigned reads_per_template = 2,
    bool strip_previous = false);

std::vector<end_t> template_ends(
    const std::vector<std::vector<std::string>>& qname_group);

void update_dup_status(std::vector<std::string>& read, bool set = true);

void write(std::ostream& out, const std::vector<std::string>& sam_record);

}
#endif // STREAMMD_MARKDUPS_H_
