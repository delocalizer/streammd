#include <cmath>
#include <iostream>
#include <string>
#include "bloomfilter.h"

#include <argparse/argparse.hpp>

#include <unistd.h>

std::string random_string( size_t length )
{
    auto randchar = []() -> char
    {
        const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
    std::string str(length,0);
    std::generate_n( str.begin(), length, randchar );
    return str;
}

int main(int argc, char* argv[]) {
  // initializes all to 0
  int maxn { 1000000 };
  float p { 0.000001 };
  std::vector<std::string> randstring(maxn);
  for (int i=0; i < maxn; i++) {
    randstring[i] = random_string(10);
  }
  std::cout << "done generating random strings" << std::endl;
  bloomfilter::BloomFilter bf(maxn, p);
  std::string item { "foo bar baz" };
  std::cout << "bf & item before add: " << (bf & item) << std::endl;
  bf |= item; 
  std::cout << "bf & item after add: " << (bf & item) << std::endl;
  for (int i=0; i < maxn; i++) {
    bf |= randstring[i];
  }
  std::cout << "count_estimate: " << bf.count_estimate() << std::endl;
}

