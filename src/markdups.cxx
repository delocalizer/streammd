#include <fstream>
#include <spdlog/spdlog.h>

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
    bool strip_previous,
    bool remove_duplicates) {

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
          process_qname_group(qname_group, out, bf, n_tpl_dup, n_aln_dup,
              reads_per_template, strip_previous, remove_duplicates);
          qname_group.clear();
        }
        qname_prev = samrecord.qname();
      }
      qname_group.push_back(std::move(samrecord));
    }
  }
  // handle the last group
  process_qname_group(qname_group, out, bf, n_tpl_dup, n_aln_dup,
      reads_per_template, strip_previous, remove_duplicates);
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
    bool strip_previous,
    bool remove_duplicates) {

  // calculate ends
  auto ends = template_ends(qname_group);

  if (ends.size() != reads_per_template) {
    std::string err = (reads_per_template == 1)
      ? "Input is not single reads?"
      : "Input is not paired or not qname-grouped?";
    throw std::runtime_error(
        qname_group[0].qname() + ": got "
        + std::to_string(ends.size()) + " primary alignment(s). " + err);
  }

  std::string signature {
    (reads_per_template == 1) 
      ? ends.front()
      : ends.front() + '_' + ends.back() };

  if (ends.front() == unmapped) {
    // sort order => if 1st is unmapped all are => do nothing.
  } else if (!bf.add(signature)) {
    // ends already seen => dupe
    n_tpl_dup += 1;
    // mark all reads in the group as dupes including secondary, supplementary
    // and unmapped: this is what Picard MD does when operating on QNAME grouped
    // (but not COORD sorted) records.
    for (auto & read : qname_group) {
      n_aln_dup += 1;
      read.update_dup_status(true);
    }
    if (remove_duplicates) {
      return;
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

// Calculate ordered template ends from the primary alignments of the qname
// group. Unmapped reads sort last by construction.
std::deque<std::string> template_ends(
    const std::vector<SamRecord>& qname_group) {

  // End strings are constructed like: "{rname}{F|R}{pos}" with the F|R between
  // because rname can end in a digit and we need string values that distinguish
  // between ("chr1", 1234) and ("chr11", 234).
  std::deque<std::string> ends;
  std::string rname, rname_prev { unmapped };
  int32_t pos, pos_prev { posmax };
  char orient;

  // Profiling shows most of the time we're doing bitarray ops (not hashing,
  // which seems surprising) but if there's more performance to be had elsewhere 
  // this is where to find it — we're not particularly clever about how we make
  // the template ends signature.
  for (auto read : qname_group) {
    // use only primary alignments for end calculation
    if ((read.flag() & flag_secondary) || (read.flag() & flag_supplementary)){
      continue;
    // unmapped
    } else if (read.flag() & flag_unmapped){
      ends.emplace_back(unmapped);
      rname_prev = unmapped;
      pos_prev = posmax;
    } else {
      rname = read.rname();
      if (read.flag() & flag_reverse){
        orient = 'R';
        pos = read.end_pos();
      } else {
        orient = 'F';
        pos = read.start_pos();
      }
      if ((rname > rname_prev) ||
          (rname == rname_prev && pos >= pos_prev)) {
        // This suggests fmt is faster than std::to_string:
        // https://www.zverovich.net/2013/09/07/integer-to-string-conversion-in-cplusplus.html
        // and that's borne out by gprof
        ends.emplace_back(rname + orient + fmt::format_int(pos).c_str());
      } else {
        ends.emplace_front(rname + orient + fmt::format_int(pos).c_str());
      }
      rname_prev = rname;
      pos_prev = pos;
    }
  }
  return ends;
}

// Log and write metrics to file as JSON.
void write_metrics(std::string metricsfname, metrics metrics) {

  float dup_frac {float(metrics.templates_marked_duplicate)/metrics.templates}; 
  spdlog::info("alignments seen: {}", metrics.alignments);
  spdlog::info("alignments marked duplicate: {}", metrics.alignments_marked_duplicate);
  spdlog::info("templates seen: {}", metrics.templates);
  spdlog::info("templates marked duplicate: {}", metrics.templates_marked_duplicate);
  spdlog::info("template duplicate fraction: {:.4f}", dup_frac);

  std::ofstream metricsf;
  metricsf.open(metricsfname);
  if (!metricsf) {
    throw std::runtime_error(metricsfname + " cannot be opened for writing");
  }
  metricsf << fmt::format(
    "{{"
      "\"ALIGNMENTS\":{},"
      "\"ALIGNMENTS_MARKED_DUPLICATE\":{},"
      "\"TEMPLATES\":{},"
      "\"TEMPLATES_MARKED_DUPLICATE\":{},"
      "\"TEMPLATE_DUPLICATE_FRACTION\":{:.4f}"
    "}}", 
    metrics.alignments,
    metrics.alignments_marked_duplicate,
    metrics.templates,
    metrics.templates_marked_duplicate,
    dup_frac);
}

}
