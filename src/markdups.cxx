#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <argparse/argparse.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/cfg/env.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "bloomfilter.h"
#include "markdups.h"
#include "version.h"


using namespace markdups;

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

// Mark duplicates in-place.
void mark_duplicates(
    std::vector<std::vector<std::string>>& qname_group,
    bloomfilter::BloomFilter& bf) {
}

// Write out records.
void write(std::vector<std::vector<std::string>>& qname_group,
    std::ostream& out) {
  for (auto record: qname_group) {
    std::string line;
    join(record, SAM_delimiter, line);
    out << line << std::endl;
  }
}

// Orchestrate the reading, duplicate marking and writing.
void process(
    std::istream& in,
    std::ostream& out,
    int reads_per_template,
    bloomfilter::BloomFilter& bf) {
  
  uint64_t n_qname { 0 };
  std::string qname_prev { "" };
  std::string qname;
  std::vector<std::vector<std::string>> qname_group;

  for (std::string line; std::getline(in, line);) {
    if (line.rfind("@", 0) == 0) {
      out << line << std::endl;
    } else {
      std::vector<std::string> fields;
      split(line, fields);
      qname = fields[0];
      if (qname != qname_prev) {
        n_qname += 1;
        if (n_qname % log_interval == 0) {
          spdlog::debug("qnames read: {0}", n_qname); 
        }
        if (!qname_group.empty()){
          mark_duplicates(qname_group, bf);
          write(qname_group, out);
        }
        qname_group.clear();
        qname_prev = qname;
      }
      qname_group.push_back(fields);
    }
  }
  mark_duplicates(qname_group, bf);
  write(qname_group, out);
}

int main(int argc, char* argv[]) {

  // https://spdlog.docsforge.com/v1.x/0.faq/#switch-the-default-logger-to-stderr
  spdlog::set_default_logger(spdlog::stderr_color_st("tmp"));
  spdlog::set_default_logger(spdlog::stderr_color_st(""));
  spdlog::cfg::load_env_levels();

  argparse::ArgumentParser cli("streammd", STREAMMD_VERSION);

  cli.add_description(
    "Read a SAM file from STDIN, mark duplicates in a single pass and stream "
    "processed records to STDOUT. Input must begin with a valid SAM header "
    "followed by qname-grouped records. Default log level is 'info' — set to "
    "something else (e.g. 'debug') via SPDLOG_LEVEL environment variable.");

  cli.add_argument("--input")
    .help("Input file. [default: STDIN]")
    .metavar("INPUT");

  cli.add_argument("--output")
    .help("Output file. [default: STDOUT]")
    .metavar("OUTPUT");

  cli.add_argument("-n", "--n-items")
    .help("Expected maximum number of templates n.")
    .default_value(default_n)
    .metavar("N_ITEMS")
    .scan<'d', uint64_t>();

  cli.add_argument("-p", "--fp-rate")
    .help("Target maximum false positive rate when n items are stored.")
    .default_value(default_p)
    .metavar("FP_RATE")
    .scan<'g', float>();

  cli.add_argument("--metrics")
    .help("Output metrics file.")
    .default_value(default_metrics)
    .metavar("METRICS_FILE");

  cli.add_argument("--single")
    .help("Accept single-ended reads as input. [default: paired-end]")
    .default_value(false)
    .implicit_value(true);

  cli.add_argument("--strip-previous")
    .help("Unset duplicate flag for any reads that have it set and are no "
          "longer considered duplicate. Only ever required if records have "
          "previously been through a duplicate marking step. [default: false]")
    .default_value(false)
    .implicit_value(true);

  try {
    cli.parse_args(argc, argv);
  } catch(const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << cli;
    std::exit(1);
  }

  auto n { cli.get<uint64_t>("-n") };
  auto p { cli.get<float>("-p") };
  auto bf { bloomfilter::BloomFilter(n, p) };

  auto infname = cli.present("--input");
  auto outfname = cli.present("--output");
  std::ifstream inf;
  std::ofstream outf;

  process(
      infname ? [&]() -> std::istream& {
        inf.open(*infname);
        if (!inf) {
          std::cerr << *infname << ": No such file" << std::endl;
          exit(1);
        }
        return inf;
      }() : std::cin,
      outfname ? [&]() -> std::ostream& { outf.open(*outfname); return outf; }() : std::cout,
      cli.get<bool>("--single") ? 1 : 2,
      bf
  );
}
