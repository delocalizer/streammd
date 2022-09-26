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

void readends(
    std::vector<std::vector<std::string>>& qname_group,
    std::vector<std::tuple<std::string, uint32_t, char>>& ends) {

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
        int num { stoi((*iter)[1]) };
        const char op { std::string((*iter)[2]).c_str()[0] };
        if (consumes_reference.count(op)) {
          ref_end += num;
        }
        ++iter;
      }
      ends.push_back(std::make_tuple(rname, ref_end + trailing_s, 'R'));
    }
  }
  sort(ends.begin(), ends.end());
}

// Mark duplicates in-place.
void mark_duplicates(
    std::vector<std::vector<std::string>>& qname_group,
    int& reads_per_template,
    bloomfilter::BloomFilter& bf){
  std::vector<std::tuple<std::string, uint32_t, char>> ends;
  readends(qname_group, ends);
  for (auto end : ends) {
    std::cout << std::get<0>(end) << "_" << std::get<1>(end) << "_" << std::get<2>(end) << std::endl;
  }
}

// Write out records.
void write(
  std::vector<std::vector<std::string>>& qname_group,
    std::ostream& out) {
  for (auto record : qname_group) {
    std::string line;
    join(record, SAM_delimiter, line);
    out << line << std::endl;
  }
}

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
      std::vector<std::string> fields;
      split(line, fields);
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

  std::string teststr1 {"10S65M"};
  std::string teststr2 {"65M10S"};
  std::smatch sm;
  regex_search(teststr1, sm, re_leading_s);
  if (!sm.empty()) {
    std::cout << sm[0] << std::endl;
  } else {
    std::cout << "not found in " << teststr1 << std::endl;
  }
  regex_search(teststr2, sm, re_trailing_s);
  if (!sm.empty()) {
    std::cout << sm[0] << std::endl;
  } else {
    std::cout << "not found in " << teststr2 << std::endl;
  }
  std::sregex_iterator iter(teststr1.begin(), teststr1.end(), re_cigar);
  std::sregex_iterator end;
  while (iter != end){
    std::cout << (*iter)[0] << std::endl;
    for (auto i = 1; i <= 2; ++i) {
      std::cout << "\t" << (*iter)[i] << std::endl;
    }
    ++iter;
  }

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
