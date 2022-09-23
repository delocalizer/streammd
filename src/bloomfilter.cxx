#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <math.h>

int main(int argc, char* argv[]) {
    // initializes all to 0
    boost::dynamic_bitset<> B1(pow(1024, 2));
    std::cout << "test(0) before set: " << B1.test(0) << std::endl;
    B1.set(0);
    std::cout << "test(0) after set: " << B1.test(0) << std::endl;
}
