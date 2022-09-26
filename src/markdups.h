#ifndef STREAMMD_MARKDUPS_H_
#define STREAMMD_MARKDUPS_H_

#include <cstdint>
#include <regex>
#include <string>
#include <tuple>

namespace markdups {
  // DEL sorts last in ASCII
  const char DEL { 127 };
  const char SAM_delimiter { '\t' };
  const float default_p { 0.000001 };
  const std::regex re_cigar { R"((?:(\d+)([MIDNSHPX=])))" };
  const std::regex re_leading_s { R"(^(\d+)S)" };                       
  const std::regex re_trailing_s { R"((\d+)S$)" };
  const std::string default_metrics { "streammd-metrics.json" };
  const std::tuple<std::string, uint32_t, char> unmapped { std::string(1, DEL), 0, DEL };
  const uint32_t log_interval { 1000000 };
  const uint64_t default_n { 1000000000 };
  const unsigned short flag_unmapped = 4;
  const unsigned short flag_reverse = 16;
  const unsigned short flag_secondary = 256;
  const unsigned short flag_duplicate = 1024;
  const unsigned short flag_supplementary = 2048;
  inline void join(
      const std::vector<std::string> fields, const char& sep, std::string& result) {
    unsigned i {0};
    auto imax { fields.size()-1 };
    for (auto & field : fields) {
      result += field;
      if (i < imax) {
        result += sep;
      }
      i++;
    }
  }

  inline void split(const std::string& line, std::vector<std::string>& fields) {
    std::stringstream line_stream{ line };
    std::string tkn;
    while(std::getline(line_stream, tkn, SAM_delimiter)) {
      fields.push_back(tkn);
    }
  }
}

#endif // STREAMMD_MARKDUPS_H_
