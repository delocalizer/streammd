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

SamRecord::SamRecord() {
  std::string buffer;
  buffer.reserve(1024);
}

// Read what we need and no more.
void SamRecord::parse() {
  if (buffer.empty()) {
    spdlog::warn("Cannot parse empty buffer");
    return;
  }
  // Seatbelts unfastened and we're just assuming good SAM records...
  size_t start, stop;

  // qname
  start=0; stop=buffer.find(SAM_delimiter, start);
  qname_ = buffer.substr(start, stop - start);

  // flag
  start=stop+1; stop=buffer.find(SAM_delimiter, start);
  flagidx_ = start;                         // needed for update
  flaglen_ = stop - start;                  // needed for update
  flag_ = stoi(buffer.substr(flagidx_, flaglen_));

  // rname
  start=stop+1; stop=buffer.find(SAM_delimiter, start);
  rname_ = buffer.substr(start, stop - start);

  // pos
  start=stop+1; stop=buffer.find(SAM_delimiter, start);
  pos_ = stoi(buffer.substr(start, stop - start));

  // cigar
  start=stop+1; stop=buffer.find(SAM_delimiter, start);
  start=stop+1; stop=buffer.find(SAM_delimiter, start);
  cigaridx_ = start;
  cigarlen_ = stop - start;
  cigar_ = buffer.substr(start, stop - start);

  // pg
  start=stop+1; stop=buffer.find(pgtag_, start);
  if (stop == std::string::npos) {
    pgidx_ = buffer.length();
    pglen_ = 0;
  } else {
    pgidx_ = stop;                          // needed for update
    stop=buffer.find(SAM_delimiter, pgidx_ + 1);
    pglen_ = (stop == std::string::npos)    // needed for update
             ? buffer.length() - pgidx_
             : stop - pgidx_;
  }
}

// Update the duplicate flag and PG:Z tag value on a record.
// When 'set' is true the flag is set, otherwise it is unset.
void SamRecord::update_dup_status(bool set) {
  auto prev = flag_;
  flag_ = set
    ? (flag_ | flag_duplicate)
    : (flag_ ^ flag_duplicate);
  if (flag_ != prev) {
    // update the buffer from the end (PG first then FLAG) so the parsed idx
    // values are valid
    buffer.replace(pgidx_, pglen_, pgtagval_);
    buffer.replace(flagidx_, flaglen_, std::to_string(flag_));
  }
}

// Return read end as tuple of the form (rname, (start|end), [FR])
end_t SamRecord::read_end() {
  // use only primary alignments for end calculation
  if ((flag_ & flag_secondary) || (flag_ & flag_supplementary)){
    return empty;
  }
  // unmapped
  if (flag_ & flag_unmapped){
    return unmapped;
  }
  // forward
  if (!(flag_ & flag_reverse)) {
    //std::smatch sm;
    //regex_search(cigar_, sm, re_leading_s);
    //int leading_s = sm.empty() ? 0 : stoi(sm[1]);
    //return std::make_tuple(rname_, pos_ - leading_s, 'F');
    return std::make_tuple(rname_, start_pos(), 'F');
  // reverse
  } else {
    std::smatch sm;
    regex_search(cigar_, sm, re_trailing_s);
    int trailing_s = sm.empty() ? 0 : stoi(sm[1]);
    int ref_end { pos_ };
    std::sregex_iterator iter(cigar_.begin(), cigar_.end(), re_cigar);
    std::sregex_iterator end;
    while (iter != end){
      if (consumes_reference.count((*iter)[2])) {
        ref_end += stoi((*iter)[1]);
      }
      ++iter;
    }
    return std::make_tuple(rname_, ref_end + trailing_s, 'R');
  }
}

// Return reference start of a fwd read (pos - leading soft clips)
int32_t SamRecord::start_pos(){
  int num { 0 };
  char ch { '\0' }, op { '\0' };
  for (size_t i=cigaridx_; i < cigaridx_ + cigarlen_; ++i) {
    ch = buffer[i] - '0';
    if (0 <= ch && ch <= 9) {
      num = 10 * num + ch;
    } else  {
      op = ch;
      break;
    }
  }
  return pos_ - ((op == 'S') ? num : 0);
}

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
      out << samrecord.buffer + "\n";
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

// Process a qname group of records; each record a vector of string fields.
void process_qname_group(
    std::vector<SamRecord>& qname_group,
    std::ostream& out,
    bloomfilter::BloomFilter& bf,
    size_t reads_per_template,
    bool strip_previous) {

  // calculate ends                               // 27 seconds [old]
  std::vector<end_t> ends { template_ends(qname_group) };

  if (ends.size() != reads_per_template) {
    std::string err = (reads_per_template == 1)
      ? "{0}: expected 1 primary alignment, got {1}. Input is not single reads?"
      : "{0}: expected 2 primary alignments, got {1}. Input is not paired or not qname-grouped?";
    spdlog::error(err, qname_group[0].qname(), ends.size());
    exit(1);
  }

  std::string ends_str { ends_to_string(ends) };
  if (ends[0] == unmapped) {
    // sort order => if 1st is unmapped all are, so do nothing.
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
    record.buffer.append(1, '\n');
    out << record.buffer;
  }
}

// Calculate template ends from the primary alignments of the qname group.
// Ends are returned as tuples of the form (rname, (start|end), [FR]) and the
// vector of them is in coordinate-sorted order. Unmapped ends sort last by
// construction.
std::vector<end_t> template_ends(
    const std::vector<SamRecord>& qname_group) {
  std::vector<end_t> ends;
  for (auto read : qname_group) {
    end_t end {read.read_end()};
    if (end == empty) {
      continue;
    }
    if (!ends.size() ||
        (std::get<0>(end) >= std::get<0>(ends[0]) &&
         std::get<1>(end) >= std::get<1>(ends[0]))) {
      ends.push_back(std::move(end));
    } else {
      ends.insert(ends.begin(), std::move(end));
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
