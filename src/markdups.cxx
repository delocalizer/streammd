#include <stdio.h>
#include <stdlib.h>

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

using end_t = std::tuple<std::string, uint32_t, char>;

namespace markdups {

// Process the input stream. Header lines are written directly to the output
// stream; reads are dispatched in qname groups for further processing. 
void process_input_stream(
    std::istream& in,
    std::ostream& out,
    bloomfilter::BloomFilter& bf,
    std::vector<std::string> cli_args,
    unsigned reads_per_template,
    bool strip_previous) {
  
  bool header_done { false };
  uint64_t n_qname { 0 };
  std::string qname_prev { "" };
  std::string qname;
  std::string header_last { "" };
  std::vector<std::vector<std::string>> qname_group;

  // Use C stdio.h getline to fill a char buffer and output with C fputs
  /* this takes ~ 5 seconds reading only; ~ 22 seconds r+w;
  char *buffer;
  size_t bufsize = 4096;
  int characters;
  buffer = (char *)malloc(bufsize * sizeof(char));
  while ((characters = getline(&buffer, &bufsize, stdin)) > 0) {
    // note c getline grabs \n also so we don't add it here
    // fputs(buffer, stdout);
  }
  */

  // Use C++ istream::getline to fill a char buffer and output with C fputs
  /* this also takes ~ 5 seconds reading only; ~ 22 seconds r+w.
  int maxbufsz {4096};
  for (char buf[maxbufsz]; in.getline(buf, maxbufsz);) {
    fputs(buf, stdout);
    fputs("\n", stdout);
  }
  if (in.eof()) {
    spdlog::trace("EOF");
  } else {
    spdlog::error("Malformed input — line too long?");
  }
  */

  // Use C++ istream::getline to fill a char buffer and output with C++ <<
  /* This took ~ 5 seconds reading only but 72(!) seconds r+w.
  int maxbufsz {4096};
  for (char buf[maxbufsz]; in.getline(buf, maxbufsz);) {
    out << buf << std::endl;
  }
  if (in.eof()) {
    spdlog::trace("EOF");
  } else {
    spdlog::error("Malformed input — line too long?");
  }
  */

  // Use C++ istream::getline to create std::string and output with C++ << 
  /* This took ~ 5 seconds reading only but 72(!) seconds r+w.
  for (std::string line; std::getline(in, line);) {
    out << line << std::endl;
  }
  */

  for (std::string line; std::getline(in, line);) {
    if (line.rfind("@", 0) == 0) {
      out << line << std::endl;
      header_last = line;
    } else {
      if (!header_done) {
        pgline(out, header_last, cli_args);
        header_done = true;
      }
      std::vector<std::string> fields { split(line, SAM_delimiter) };
      qname = fields[0];
      if (qname != qname_prev) {
        n_qname += 1;
        if (n_qname % log_interval == 0) {
          spdlog::debug("qnames read: {0}", n_qname); 
        }
        if (qname_group.size()) {
          process_qname_group(qname_group, out, bf, reads_per_template, strip_previous);
          qname_group.clear();
        }
        qname_prev = qname;
      }
      qname_group.emplace_back(fields);
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
  std::string cl { "CL:" + join(cli_args, ' ') };
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
  out << join(tags, '\t') << std::endl;
}

// Process a qname group of records; each record a vector of string fields.
void process_qname_group(
    std::vector<std::vector<std::string>>& qname_group,
    std::ostream& out,
    bloomfilter::BloomFilter& bf,
    unsigned reads_per_template,
    bool strip_previous) {

  // calculate ends
  std::vector<end_t> ends { template_ends(qname_group) };

  if (ends.size() != reads_per_template) {
    std::string err = (reads_per_template == 1)
      ? "{0}: expected 1 primary alignment, got {1}. Input is not single reads?"
      : "{0}: expected 2 primary alignments, got {1}. Input is not paired or not qname-grouped?";
    spdlog::error(err, qname_group[0][0], ends.size());
    exit(1);
  }

  std::string ends_str { ends_to_string(ends) };
  if (ends[0] == unmapped) {
    // sort order => if 1st is unmapped all are, so do nothing.
  } else if (!bf.add(ends_str)) {
    // ends already seen => dupe
    for (auto & read : qname_group) {
      update_dup_status(read, true);
    }
  } else if (strip_previous) {
    for (auto & read : qname_group) {
      update_dup_status(read, false);
    }
  }
  // write to output
  for (auto record : qname_group) {
    out << join(record, SAM_delimiter) << std::endl;
  }
}

// Calculate template ends from the primary alignments of the qname group.
// Ends are returned as tuples of the form (rname, (start|end), [FR]) and the
// vector of them is in coordinate-sorted order. Unmapped ends sort last by
// construction.
std::vector<end_t> template_ends(
    const std::vector<std::vector<std::string>>& qname_group) {

  std::vector<end_t> ends;

  for (auto read : qname_group) {
    auto flag      = stoi(read[1]);
    auto rname     = read[2];
    auto ref_start = stoi(read[3]);
    auto cigar     = read[5];
    // use only primary alignments for end calculation
    if ((flag & flag_secondary) || (flag & flag_supplementary)){
      continue;
    }
    // unmapped
    if (flag & flag_unmapped){
      ends.emplace_back(unmapped);
    // forward
    } else if (!(flag & flag_reverse)) {
      std::smatch sm;
      regex_search(cigar, sm, re_leading_s);
      int leading_s = sm.empty() ? 0 : stoi(sm[1]);
      ends.emplace_back(std::make_tuple(rname, ref_start - leading_s, 'F'));
    // reverse
    } else {
      std::smatch sm;
      regex_search(cigar, sm, re_trailing_s);
      int trailing_s = sm.empty() ? 0 : stoi(sm[1]);
      int ref_end { ref_start };
      std::sregex_iterator iter(cigar.begin(), cigar.end(), re_cigar);
      std::sregex_iterator end;
      while (iter != end){
        if (consumes_reference.count((*iter)[2])) {
          ref_end += stoi((*iter)[1]);
        }
        ++iter;
      }
      ends.emplace_back(std::make_tuple(rname, ref_end + trailing_s, 'R'));
    }
  }
  sort(ends.begin(), ends.end());
  return ends;
}

// Update the duplicate flag and PG:Z tag value on a read.
// When 'set' is true the flag is set, otherwise it is unset.
void update_dup_status(std::vector<std::string>& read, bool set) {
  auto prev = read[1];
  auto flag = stoi(prev);
  // update the flag
  read[1] = set
    ? std::to_string(flag | flag_duplicate)
    : std::to_string(flag ^ flag_duplicate);
  // update PG:Z if flag changed
  if (read[1] != prev) {
    auto imax { read.size() - 1 };
    std::string pg_old;
    for (auto i = sam_opts_idx; i <= imax; i++) {
      if (read[i].rfind(pgtag, 0) == 0) {
        pg_old = read[i];
        read[i] = pgtag_val;
      }
    }
    if (pg_old.empty()) {
      read.emplace_back(pgtag_val);
    }
  }
}

}

using namespace markdups;

int main(int argc, char* argv[]) {

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
  std::ios::sync_with_stdio(false);

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
