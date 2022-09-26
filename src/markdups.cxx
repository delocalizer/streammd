#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include <argparse/argparse.hpp>

#include "bloomfilter.h"
#include "markdups.h"
#include "version.h"


// Set logging filter from LOG_LEVEL environment variable, defaulting to info.
void logging_init()
{
  namespace logging = boost::log;
  static std::map<std::string, logging::trivial::severity_level> const levels {
    {
      {"trace",   logging::trivial::severity_level::trace},
      {"debug",   logging::trivial::severity_level::debug},
      {"info",    logging::trivial::severity_level::info},
      {"warning", logging::trivial::severity_level::warning},
      {"error",   logging::trivial::severity_level::error},
      {"fatal",   logging::trivial::severity_level::fatal}
    }
  };
  std::string severity_str { std::getenv("LOG_LEVEL") ?: "info" };
  logging::trivial::severity_level severity;
  auto entry = levels.find(severity_str);
  if (entry != levels.end()) {
    severity = entry->second;
  } else {
    std::cerr << severity_str << " is not a valid log level" << std::endl;
    std::exit(1);
  }
  logging::core::get()->set_filter
  (
    logging::trivial::severity >= severity
  );
}

using namespace markdups;

void read_input(std::istream& in, std::ostream& out) {
  for (std::string line; std::getline(in, line);) {
    if (line.rfind("@", 0) == 0) {
      out << line << std::endl;
    }
  }
}

int main(int argc, char* argv[]) {

  logging_init();

  argparse::ArgumentParser cli("streammd", STREAMMD_VERSION);

  cli.add_description(
    "Read a SAM file from STDIN, mark duplicates in a single pass and stream "
    "processed records to STDOUT. Input must begin with a valid SAM header "
    "followed by qname-grouped records. Default log level is 'info' — set to "
    "something else (e.g. 'warning') with the LOG_LEVEL environment variable.");

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
  auto reads_per_template { cli.get<bool>("--single") ? 1 : 2 };
  auto infname = cli.present("--input");
  auto outfname = cli.present("--output");
  std::ifstream inf;
  std::ofstream outf;
  read_input(
      infname ? [&]() -> std::istream& {
        inf.open(*infname);
        if (!inf) {
          std::cerr << *infname << " No such file" << std::endl;
          exit(1);
        }
        return inf;
      }() : std::cin,
      outfname ? [&]() -> std::ostream& { outf.open(*outfname); return outf; }() : std::cout
  );

  //auto bf { bloomfilter::BloomFilter(n, p) };
}

