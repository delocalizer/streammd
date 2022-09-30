#ifndef STREAMMD_MARKDUPS_H_
#define STREAMMD_MARKDUPS_H_

#include <regex>

#include "bloomfilter.h"

namespace markdups {

// DEL sorts last in ASCII
const char DEL { 127 };
const char SAM_delimiter { '\t' };
const float default_p { 0.000001 };
const std::regex re_pgid { R"(\tID:([^\t]+))" };
const std::string pgid { "streammd"  };
const std::string pgtag { "PG:Z:" };
const std::string default_metrics { pgid + "-metrics.json" };
const std::string unmapped { DEL };
const uint32_t log_interval { 1000000 };
const uint64_t default_n { 1000000000 };
const size_t short flag_unmapped      = 4;
const size_t short flag_reverse       = 16;
const size_t short flag_secondary     = 256;
const size_t short flag_duplicate     = 1024;
const size_t short flag_supplementary = 2048;

// A highly specialized representation of a SAM record, for the purposes of
// duplicate marking — probably not generally useful for other things...
class SamRecord {

 public:
  SamRecord() {
    std::string buffer;
    buffer.reserve(1024);
  }
  std::string buffer;
  inline const std::string& qname() { return qname_; };
  inline const uint16_t& flag() { return flag_; };
  inline const std::string& rname() { return rname_; };

  // Construct what we need and no more.
  inline void parse() {
    if (buffer.empty()) {
      spdlog::warn("Cannot parse empty buffer");
      return;
    }
    size_t start, stop;

    // Seatbelts unfastened and we're just assuming good SAM records...

    // qname
    start=0; stop=buffer.find(SAM_delimiter, start);
    qname_ = buffer.substr(start, stop - start);

    // flag
    start=stop+1; stop=buffer.find(SAM_delimiter, start);
    flagidx_ = start;
    flaglen_ = stop - start;
    flag_ = stoi(buffer.substr(flagidx_, flaglen_));

    // rname
    start=stop+1; stop=buffer.find(SAM_delimiter, start);
    rname_ = buffer.substr(start, stop - start);

    // pos
    start=stop+1; stop=buffer.find(SAM_delimiter, start);
    pos_ = stoi(buffer.substr(start, stop - start));

    // cigar
    start=stop+1; stop=buffer.find(SAM_delimiter, start);
    start=stop+1; stop=buffer.find(SAM_delimiter, start);
    cigaridx_ = start;
    cigarlen_ = stop - start;

    // pg
    start=stop+1; stop=buffer.find(pgtag_, start);
    if (stop == std::string::npos) {
      pgidx_ = buffer.length();
      pglen_ = 0;
    } else {
      pgidx_ = stop;
      stop=buffer.find(SAM_delimiter, pgidx_ + 1);
      pglen_ = (stop == std::string::npos)
               ? buffer.length() - pgidx_
               : stop - pgidx_;
    }
  }

  // Return reference start of a fwd read (pos - leading soft clips)
  inline int32_t start_pos(){
    int num { 0 };
    char ch { '\0' }, op { '\0' };
    for (size_t i=cigaridx_; i < cigaridx_ + cigarlen_; ++i) {
      ch = buffer[i] - '0';
      if (0 <= ch && ch <= 9) {
        num = 10 * num + ch;
      } else  {
        op = buffer[i];
        break;
      }
    }
    return pos_ - ((op == 'S') ? num : 0);
  }

  // Return reference end of a rev read (pos + reflen + trailing soft clips)
  inline int32_t end_pos(){
    int num { 0 }, prev { 0 }, reflen { 0 };
    char ch { '\0' }, op { '\0' };
    for (size_t i=cigaridx_; i < cigaridx_ + cigarlen_; ++i) {
      ch = buffer[i] - '0';
      if (0 <= ch && ch <= 9) {
        num = 10 * num + ch;
      } else  {
        op = buffer[i];
        // ops that consume reference
        if (op == 'M' ||
            op == 'D' ||
            op == 'N' ||
            op == '=' ||
            op == 'X') {
          reflen += num;
        }
        prev = num;
        num = 0;
      }
    }
    return pos_ + reflen + ((op == 'S') ? prev : 0);
  }

  // Update the duplicate flag and PG:Z tag value on a record.
  // When 'set' is true the flag is set, otherwise it is unset.
  inline void update_dup_status(bool set) {
    auto prev = flag_;
    flag_ = set
      ? (flag_ | flag_duplicate)
      : (flag_ ^ flag_duplicate);
    if (flag_ != prev) {
      // update the buffer from the end (PG first then FLAG) so the parsed idx
      // values are valid
      buffer.replace(pgidx_, pglen_, pgtagval_);
      buffer.replace(flagidx_, flaglen_, std::to_string(flag_));
    }
  }

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
  size_t pgidx_;
  size_t pglen_;
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

std::deque<std::string> template_ends(
    const std::vector<SamRecord>& qname_group);

}
#endif // STREAMMD_MARKDUPS_H_
