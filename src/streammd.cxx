#include <fstream>
#include <iostream>

#include <argparse/argparse.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/cfg/env.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "bloomfilter.h"
#include "markdups.h"
#include "version.h"

using namespace markdups;

int main(int argc, char* argv[]) {

  // Un-sync and un-tie the I/O; this makes a HUGE difference in speed when
  // reading from stdin and writing to stdout with streams.
  std::ios::sync_with_stdio(false);
  std::cin.tie(0);

  spdlog::set_default_logger(spdlog::stderr_color_st("main"));
  spdlog::cfg::load_env_levels();

  argparse::ArgumentParser cli(pgid, STREAMMD_VERSION);

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
    .scan<'g', double>();

  cli.add_argument("--allow-overcapacity")
    .help("Warn instead of error when Bloom filter capacity is exceeded. "
          "[default: false]")
    .default_value(false)
    .implicit_value(true);

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
    spdlog::error(err.what());
    std::exit(1);
  }

  auto n { cli.get<uint64_t>("-n") };
  auto p { cli.get<double>("-p") };
  auto bf { bloomfilter::BloomFilter(p, n) };
  std::vector<std::string> args(argv, argv + argc);

  auto infname = cli.present("--input");
  auto outfname = cli.present("--output");
  auto metricsfname = cli.get("--metrics");
  std::ifstream inf;
  std::ofstream outf;

  try {
    auto result = process_input_stream(
        infname ? [&]() -> std::istream& {
          inf.open(*infname);
          if (!inf) {
            std::string errmsg { *infname + " cannot be opened for reading" };
            throw std::runtime_error(errmsg);
          }
          return inf;
        }() : std::cin,
        outfname ? [&]() -> std::ostream& {
          outf.open(*outfname);
          if (!outf) {
            std::string errmsg { *outfname + " cannot be opened for writing" };
            throw std::runtime_error(errmsg);
          }
          return outf;
        }() : std::cout,
        bf,
        args,
        cli.get<bool>("--single") ? 1 : 2,
        cli.get<bool>("--strip-previous")
    );
    auto cap { float(result.templates)/bf.n() };
    if (cap <= 1.0) {
      spdlog::info("BloomFilter: {:.3g}% capacity", 100*cap);
      write_metrics(metricsfname, result);
    } else if (cli.get<bool>("--allow-overcapacity")) {
      write_metrics(metricsfname, result);
      spdlog::warn("BloomFilter: {:.3g}% capacity", 100*cap);
      spdlog::warn("BloomFilter capacity of {} exceeded: "
        "false positive rate target of {} is likely violated.", bf.n(), bf.p());
    } else {
      spdlog::error("BloomFilter: {:.3g}% capacity", 100*cap);
      throw std::runtime_error("BloomFilter capacity exceeded");
    }
  } catch(const std::runtime_error& err) {
    spdlog::error(err.what());
    std::exit(1);
  }
}
