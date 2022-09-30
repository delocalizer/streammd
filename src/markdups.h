#ifndef STREAMMD_MARKDUPS_H_
#define STREAMMD_MARKDUPS_H_

#include <cstdint>
#include <regex>
#include <set>
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
const std::regex re_pgid { R"(\tID:([^\t]+))" };
const std::regex re_trailing_s { R"((\d+)S$)" };
const std::set<std::string> consumes_reference { "M", "D", "N", "=", "X" };
const std::string pgid { "streammd"  };
const std::string pgtag { "PG:Z:" };
const std::string default_metrics { pgid + "-metrics.json" };
const end_t empty { "", -1, '\0' };
const end_t unmapped { std::string(1, DEL), -1, DEL };
const uint32_t log_interval { 1000000 };
const uint64_t default_n { 1000000000 };
const size_t short flag_unmapped      = 4;
const size_t short flag_reverse       = 16;
const size_t short flag_secondary     = 256;
const size_t short flag_duplicate     = 1024;
const size_t short flag_supplementary = 2048;
const size_t short sam_opts_idx = 11;

// This is a very specialized representation of a SAM record, for the
// purposes of duplicate marking — not generally useful for much else
class SamRecord {

 public:
  SamRecord();
  std::string buffer;
  void parse();
  void update_dup_status(bool set);
  end_t read_end();
  inline const std::string& qname() { return qname_; };
  inline const uint16_t& flag() { return flag_; };
  inline const std::string& rname() { return rname_; };
  inline const int32_t& pos() { return pos_; };
  inline const std::string& cigar() { return cigar_; };

 private:
  inline static const std::string pgtag_ { std::string(1, SAM_delimiter) + pgtag };
  inline static const std::string pgtagval_ { pgtag_ + pgid };
  std::string qname_;
  size_t flagidx_;
  size_t flaglen_;
  uint16_t flag_;
  std::string rname_;
  int32_t pos_;
  size_t cigaridx_;
  size_t cigarlen_;
  std::string cigar_;
  size_t pgidx_;
  size_t pglen_;
  int32_t start_pos();
};

inline std::string join(
    const std::vector<std::string>& elems,
    const char sep,
    const char end = '\0') {
  std::string joined;
  auto stop { elems.end() };
  for (auto i = elems.begin(); i != stop; ++i) {
    joined += *i;
    if (i != stop - 1) { joined += sep; }
  }
  if (end != '\0') {
    joined += end;
  }
  return joined;
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
    std::vector<std::string> cli_args,
    size_t reads_per_template = 2,
    bool strip_previous = false);

void pgline(
    std::ostream& out,
    const std::string& header_last,
    const std::vector<std::string>& cli_args);

void process_qname_group(
    std::vector<SamRecord>& qname_group,
    std::ostream& out,
    bloomfilter::BloomFilter& bf,
    size_t reads_per_template = 2,
    bool strip_previous = false);

std::vector<end_t> template_ends(
    const std::vector<SamRecord>& qname_group);

}
#endif // STREAMMD_MARKDUPS_H_
