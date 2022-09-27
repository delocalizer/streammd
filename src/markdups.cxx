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

// Orchestrate the reading, dupe marking, and writing.
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
      std::vector<std::string> fields { split(line, SAM_delimiter) };
      qname = fields[0];
      if (qname != qname_prev) {
        n_qname += 1;
        if (n_qname % log_interval == 0) {
          spdlog::debug("qnames read: {0}", n_qname); 
        }
        mark_duplicates(qname_group, reads_per_template, bf);
        //write(qname_group, out);
        qname_group.clear();
        qname_prev = qname;
      }
      qname_group.push_back(fields);
    }
  }
  mark_duplicates(qname_group, reads_per_template, bf);
  //write(qname_group, out);
}

// Mark duplicates in-place.
void mark_duplicates(
    std::vector<std::vector<std::string>>& qname_group,
    int reads_per_template,
    bloomfilter::BloomFilter& bf){
  std::vector<end_t> ends { readends(qname_group) };
  std::vector<std::string> ends_str;
  for (auto end : ends) {
    std::string end_str { std::get<0>(end) };
    // Place F|R next because rname can end in a digit and we need string
    // values that distinguish between ("chr1", 1234) and ("chr11", 234) 
    end_str += std::get<2>(end);
    end_str += std::to_string(std::get<1>(end));
    ends_str.push_back(end_str);
    std::cout << end_str << std::endl;
  }
}

// Calculate read ends.
// Ends are returned as tuples of the form (rname, (start|end), [FR]) and the
// vector of them is in coordinate-sorted order. Unmapped ends sort last by
// construction.
std::vector<end_t> readends(
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
      ends.push_back(unmapped);
    // forward
    } else if (!(flag & flag_reverse)) {
      std::smatch sm;
      regex_search(cigar, sm, re_leading_s);
      int leading_s = sm.empty() ? 0 : stoi(sm[1]);
      ends.push_back(std::make_tuple(rname, ref_start - leading_s, 'F'));
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
      ends.push_back(std::make_tuple(rname, ref_end + trailing_s, 'R'));
    }
  }
  sort(ends.begin(), ends.end());
  return ends;
}

// Write out records.
void write(
  std::vector<std::vector<std::string>>& qname_group,
    std::ostream& out) {
  for (auto record : qname_group) {
    std::string line { join(record, SAM_delimiter) };
    out << line << std::endl;
  }
}

}

using namespace markdups;

int main(int argc, char* argv[]) {

  spdlog::set_default_logger(spdlog::stderr_color_st("main"));
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
