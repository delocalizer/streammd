#include <fstream>
#include <iostream>

#include <argparse/argparse.hpp>
#include <bytesize/bytesize.hh>
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

  cli.add_argument("-p", "--fp-rate")
    .help("Target maximum false positive rate.")
    .default_value(double(0.000001))
    .metavar("FP_RATE")
    .scan<'g', double>();

  cli.add_argument("-m", "--mem")
    .help("Memory allowance for the Bloom filter, e.g \"4GiB\". Both binary "
          "(kiB|MiB|GiB) and decimal (kB|MB|GB) formats are understood. "
          "As a result of implementation details, a value that is an exact "
          "power of 2 (512MiB, 1GiB, 2GiB etc) gives a modest processing "
          "speed advantage (~5%) over neighbouring values.")
    .default_value(std::string("4GiB"))
    .metavar("MEM");

  cli.add_argument("--allow-overcapacity")
    .help("Warn instead of error when Bloom filter capacity is exceeded. "
          "[default: false]")
    .default_value(false)
    .implicit_value(true);

  cli.add_argument("--metrics")
    .help("Output metrics file.")
    .default_value(std::string(pgid + "-metrics.json"))
    .metavar("METRICS_FILE");

  cli.add_argument("--remove-duplicates")
    .help("Omit detected duplicates from the output.")
    .default_value(false)
    .implicit_value(true);

  cli.add_argument("--show-capacity")
    .help("Do no work, just print the capacity of the Bloom filter that would "
          "be constructed with the given --fp-rate and --mem values.")
    .default_value(false)
    .implicit_value(true);

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
    std::vector<std::string> args(argv, argv + argc);
    auto p { cli.get<double>("-p") };
    auto mem { cli.get("-m") };

    if (cli.get<bool>("--show-capacity")){
      auto mem_actual { bloomfilter::BloomFilter::memspec_to_bytes(mem, false) };
      auto mem_actual_str { bytesize::bytesize(mem_actual).format() };
      std::cout << "specified fp-rate:       " << p << std::endl;
      std::cout << "specified mem:           " << mem << std::endl;
      std::cout << "allocated mem:           " << mem_actual_str << std::endl;
      std::cout << "Bloom filter capacity n: "
                << bloomfilter::BloomFilter::capacity(p, 8*mem_actual, 10)
                << std::endl;
      return 0;
    }

    auto bf { bloomfilter::BloomFilter::fromMemSpec(p, mem, false) };
    spdlog::info("BloomFilter capacity: {} items", bf.n());

    auto infname = cli.present("--input");
    auto outfname = cli.present("--output");
    auto metricsfname = cli.get("--metrics");
    std::ifstream inf;
    std::ofstream outf;

    auto result = process_input_stream(
        infname ? [&]() -> std::istream& {
          inf.open(*infname);
          if (!inf) {
            throw std::runtime_error(*infname + " cannot be opened for reading");
          }
          return inf;
        }() : std::cin,
        outfname ? [&]() -> std::ostream& {
          outf.open(*outfname);
          if (!outf) {
            throw std::runtime_error(*outfname + " cannot be opened for writing");
          }
          return outf;
        }() : std::cout,
        bf,
        args,
        cli.get<bool>("--single") ? 1 : 2,
        cli.get<bool>("--strip-previous"),
        cli.get<bool>("--remove-duplicates")
    );
    auto cap { float(result.templates)/bf.n() };
    if (cap <= 1.0) {
      spdlog::info("BloomFilter at {:.3g}% capacity", 100*cap);
      write_metrics(metricsfname, result);
    } else if (cli.get<bool>("--allow-overcapacity")) {
      write_metrics(metricsfname, result);
      spdlog::warn("{} items leaves BloomFilter at {:.3g}% capacity.",
                   result.templates, 100*cap);
      spdlog::warn("False positive rate {} exceeded.", bf.p());
    } else {
      spdlog::error("{} items leaves BloomFilter at {:.3g}% capacity.",
                    result.templates, 100*cap);
      throw std::runtime_error("BloomFilter capacity exceeded.");
    }
  } catch(const std::exception& err) {
    spdlog::error(err.what());
    return 1;
  }
}
