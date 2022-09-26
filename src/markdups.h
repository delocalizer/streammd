#ifndef STREAMMD_MARKDUPS_H_
#define STREAMMD_MARKDUPS_H_

#include <cstdint>
#include <string>
#include <tuple>

namespace markdups {
  const char SAM_delimiter { '\t' };
  const char DEL { 127 };
  const float default_p { 0.000001 };
  const std::string default_metrics { "streammd-metrics.json" };
  // DEL sorts last in ASCII
  const std::string unmapped { DEL };
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
