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
metrics process_input_stream(
    std::istream& in,
    std::ostream& out,
    bloomfilter::BloomFilter& bf,
    std::vector<std::string> cli_args,
    size_t reads_per_template,
    bool strip_previous) {

  bool header_done { false };
  uint64_t n_tpl { 0 }, n_aln {0}, n_tpl_dup {0}, n_aln_dup {0};
  std::string header_prev { "" };
  std::string qname;
  std::string qname_prev { "" };
  std::vector<SamRecord> qname_group;

  for (SamRecord samrecord; std::getline(in, samrecord.buffer);) {
    // header
    if (samrecord.buffer[0] == '@') {
      out << samrecord.buffer;
      out << "\n";
      header_prev = samrecord.buffer;
    } else {
      n_aln += 1;
      if (!header_done) {
        pgline(out, header_prev, cli_args);
        header_done = true;
      }
      samrecord.parse();
      if (samrecord.qname() != qname_prev) {
        n_tpl += 1;
        if (n_tpl % log_interval == 0) {
          spdlog::debug("qnames seen: {0}", n_tpl);
        }
        if (qname_group.size()) {
          process_qname_group(qname_group, out, bf,
              n_tpl_dup, n_aln_dup, reads_per_template, strip_previous);
          qname_group.clear();
        }
        qname_prev = samrecord.qname();
      }
      qname_group.push_back(std::move(samrecord));
    }
  }
  // handle the last group
  process_qname_group(qname_group, out, bf,
      n_tpl_dup, n_aln_dup, reads_per_template, strip_previous);
  return metrics { n_tpl, n_tpl_dup, n_aln, n_aln_dup };
}

// Write @PG line to the output stream; to be called after all existing header
// lines have been processed.
// 'header_prev' is the currently last line of the header, which may contain a
// previous @PG entry.
// 'cli_args' is a string vector of argv.
void pgline(
    std::ostream& out,
    const std::string& header_prev,
    const std::vector<std::string>& cli_args) {
  std::vector<std::string> tags {
    "@PG",
    "ID:" + pgid,
    "PN:" + pgid,
    "CL:" + join(cli_args, ' '),
    "VN:" + std::string(STREAMMD_VERSION)
  };
  if (header_prev.rfind("@PG\t", 0) == 0) {
    std::smatch sm;
    regex_search(header_prev, sm, re_pgid);
    if (sm.empty()) {
      spdlog::warn("{}: invalid @PG line lacks ID: tag", header_prev);
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
    uint64_t& n_tpl_dup,
    uint64_t& n_aln_dup,
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
    n_tpl_dup += 1;
    for (auto & read : qname_group) {
      n_aln_dup += 1;
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

// Log and write metrics to file as JSON.
void write_metrics(std::string metricsfname, metrics metrics) {

  float dup_frac {float(metrics.templates_marked_duplicate)/metrics.templates}; 
  spdlog::info("template duplicate fraction: {:.4f}", dup_frac);
  spdlog::info("templates seen: {}", metrics.templates);
  spdlog::info("templates marked duplicate: {}", metrics.templates_marked_duplicate);
  spdlog::info("alignments seen: {}", metrics.alignments);
  spdlog::info("alignments marked duplicate: {}", metrics.alignments_marked_duplicate);

  // std::format doesn't arrive till c++20 so format file output with ye olde
  // fprintf. Note that spdlog includes fmt lib but using that would tie us to
  // spdlog for logging.
  FILE* metricsf { fopen(metricsfname.c_str(), "w") };
  if (nullptr == metricsf) {
    spdlog::error("{} cannot be opened for writing.", metricsfname);
    exit(1);
  } else {
    fprintf(metricsf,
        "{"
          "\"ALIGNMENTS\":%lu,"
          "\"ALIGNMENTS_MARKED_DUPLICATE\":%lu,"
          "\"TEMPLATES\":%lu,"
          "\"TEMPLATES_MARKED_DUPLICATE\":%lu,"
          "\"TEMPLATE_DUPLICATE_FRACTION\":%0.4f"
        "}",
        metrics.alignments,
        metrics.alignments_marked_duplicate,
        metrics.templates,
        metrics.templates_marked_duplicate,
        dup_frac);
  }
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
  auto metricsfname = cli.get("--metrics");
  std::ifstream inf;
  std::ofstream outf;
  //int sz {4096}; // fstat reports good for file output on my laptop
  int sz {1024}; // fstat reports good for stdout on my laptop; also seems good for file redirect
  std::vector<char> buf;
  buf.resize(sz);
  outf.rdbuf()->pubsetbuf(&buf[0], sz);

  auto result = process_input_stream(
      infname ? [&]() -> std::istream& {
        inf.open(*infname);
        if (!inf) {
          spdlog::error("{}: cannot be opened for reading", *infname);
          exit(1);
        }
        return inf;
      }() : std::cin,
      outfname ? [&]() -> std::ostream& {
        outf.open(*outfname);
        if (!outf) {
          spdlog::error("{}: cannot be opened for writing", *outfname);
          exit(1);
        }
        return outf;
      }() : std::cout,
      bf,
      args,
      cli.get<bool>("--single") ? 1 : 2,
      cli.get<bool>("--strip-previous")
  );
  write_metrics(metricsfname, result);
}
