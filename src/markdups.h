#ifndef STREAMMD_MARKDUPS_H_
#define STREAMMD_MARKDUPS_H_

#include <cstdint>
#include <string>

namespace markdups {
  const uint64_t default_n { 1000000000 };
  const float default_p { 0.000001 };
  const std::string default_metrics { "streammd-metrics.json" };
}

#endif // STREAMMD_MARKDUPS_H_
