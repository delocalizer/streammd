#include <fstream>
#include <iostream>
#include <vector>

#include <argparse/argparse.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/cfg/env.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "bloomfilter.h"
#include "markdups.h"
#include "version.h"

namespace markdups {

// Process the input stream. Header lines are written directly to the output
// stream; reads are dispatched in qname groups for further processing.
void process_input_stream(
    std::istream& in,
    std::ostream& out,
    bloomfilter::BloomFilter& bf,
    std::vector<std::string> cli_args,
    size_t reads_per_template,
    bool strip_previous) {

  bool header_done { false };
  uint64_t n_qname { 0 };
  std::string qname_prev { "" };
  std::string qname;
  std::string header_last { "" };
  std::vector<SamRecord> qname_group;

  for (SamRecord samrecord; std::getline(in, samrecord.buffer);) {
    // header
    if (samrecord.buffer[0] == '@') {
      out << samrecord.buffer;
      out << "\n";
      header_last = samrecord.buffer;
    } else {
      if (!header_done) {
        pgline(out, header_last, cli_args);
        header_done = true;
      }
      samrecord.parse();
      if (samrecord.qname() != qname_prev) {
        n_qname += 1;
        if (n_qname % log_interval == 0) {
          spdlog::debug("qnames read: {0}", n_qname);
        }
        if (qname_group.size()) {
          process_qname_group(qname_group, out, bf, reads_per_template, strip_previous);
          qname_group.clear();
        }
        qname_prev = samrecord.qname();
      }
      qname_group.push_back(std::move(samrecord));
    }
  }
  // handle the last group
  process_qname_group(qname_group, out, bf, reads_per_template, strip_previous);
}

// Write @PG line to the output stream; to be called after all existing header
// lines have been processed.
// 'header_last' is the currently last line of the header, which may contain a
// previous @PG entry.
// 'cli_args' is a string vector of argv.
void pgline(
    std::ostream& out,
    const std::string& header_last,
    const std::vector<std::string>& cli_args) {
  std::vector<std::string> tags {
    "@PG",
    "ID:" + pgid,
    "PN:" + pgid,
    "CL:" + join(cli_args, ' '),
    "VN:" + std::string(STREAMMD_VERSION)
  };
  if (header_last.rfind("@PG\t", 0) == 0) {
    std::smatch sm;
    regex_search(header_last, sm, re_pgid);
    if (sm.empty()) {
      spdlog::warn("{}: invalid @PG line lacks ID: tag", header_last);
    } else {
      tags.emplace_back("PP:" + std::string(sm[1]));
    }
  }
  out << join(tags, SAM_delimiter, '\n');
}

// Process a qname group of records.
void process_qname_group(
    std::vector<SamRecord>& qname_group,
    std::ostream& out,
    bloomfilter::BloomFilter& bf,
    size_t reads_per_template,
    bool strip_previous) {

  // calculate ends
  auto ends = template_ends(qname_group);

  if (ends.size() != reads_per_template) {
    std::string err = (reads_per_template == 1)
      ? "{0}: expected 1 primary alignment, got {1}. Input is not single reads?"
      : "{0}: expected 2 primary alignments, got {1}. Input is not paired or not qname-grouped?";
    spdlog::error(err, qname_group[0].qname(), ends.size());
    exit(1);
  }

  std::string ends_str {
    (reads_per_template == 1) 
      ? ends.front()
      : ends.front() + '_' + ends.back() };

  if (ends.front() == unmapped) {
    // sort order => if 1st is unmapped all are => do nothing.
  } else if (!bf.add(ends_str)) {
    // ends already seen => dupe
    for (auto & read : qname_group) {
      read.update_dup_status(true);
    }
  } else if (strip_previous) {
    for (auto & read : qname_group) {
      read.update_dup_status(false);
    }
  }
  // write to output
  for (auto record : qname_group) {
    out << record.buffer;
  }
}

// Calculate template ends from the primary alignments of the qname group.
// Unmapped reads sort last by construction.
std::deque<std::string> template_ends(
    const std::vector<SamRecord>& qname_group) {

  // End strings are constructed like: "{rname}{F|R}{pos}" with the F|R between
  // because rname can end in a digit and we need string values that distinguish
  // between ("chr1", 1234) and ("chr11", 234).
  std::deque<std::string> ends;

  for (auto read : qname_group) {
    // use only primary alignments for end calculation
    if ((read.flag() & flag_secondary) || (read.flag() & flag_supplementary)){
      continue;
    // unmapped
    } else if (read.flag() & flag_unmapped){
      ends.emplace_back(unmapped);
    // reverse
    } else if (read.flag() & flag_reverse) {
      ends.emplace_back(read.rname() + 'R' + std::to_string(read.end_pos()));
    // forward
    } else {
      ends.emplace_front(read.rname() + 'F' + std::to_string(read.start_pos()));
    }
  }
  return ends;
}

}

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
  std::vector<std::string> args(argv, argv + argc);

  auto infname = cli.present("--input");
  auto outfname = cli.present("--output");
  std::ifstream inf;
  std::ofstream outf;

  process_input_stream(
      infname ? [&]() -> std::istream& {
        inf.open(*infname);
        if (!inf) {
          std::cerr << *infname << ": No such file" << std::endl;
          exit(1);
        }
        return inf;
      }() : std::cin,
      outfname ? [&]() -> std::ostream& { outf.open(*outfname); return outf; }() : std::cout,
      bf,
      args,
      cli.get<bool>("--single") ? 1 : 2,
      cli.get<bool>("--strip-previous")
  );
}
