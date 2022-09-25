#include <cmath>
#include <iostream>
#include <string>
#include "bloomfilter.h"

#include <argparse/argparse.hpp>

#include "streammd.h"
#include "version.h"

using namespace streammd;

int main(int argc, char* argv[]) {

  argparse::ArgumentParser cli("streammd", STREAMMD_VERSION);

  cli.add_argument("-n", "--n-items")
    .help("Expected maximum number of templates n.")
    .default_value(default_n)
    .scan<'d', uint64_t>();

  cli.add_argument("-p", "--fp-rate")
    .help("Target maximum false positive rate when n items are stored.")
    .default_value(default_p)
    .scan<'g', float>();

  try {
    cli.parse_args(argc, argv);
  } catch(const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << cli;
    std::exit(1);
  }

  auto n { cli.get<uint64_t>("-n") };
  auto p { cli.get<float>("-p") };
  std::cout << "n: " << n << " p: " << p << std::endl;
  bloomfilter::BloomFilter bf(n, p);
  std::string item { "foo bar baz" };
  std::cout << "bf.contains(item) before add: " << bf.contains(item) << std::endl;
  bf.add(item); 
  std::cout << "bf.contains(item) after add: " << bf.contains(item) << std::endl;
  std::cout << "count_estimate: " << bf.count_estimate() << std::endl;
}

