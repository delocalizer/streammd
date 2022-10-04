AM_CXXFLAGS = \
	      -Wall \
	      -Wextra \
	      -std=c++17 \
	      -I ../external \
	      -I ../src 

TESTS = driver.sh 

check_PROGRAMS = streammd_tests
streammd_tests_SOURCES = test_bloomfilter.cxx test_markdups.cxx
streammd_tests_LDADD = ../src/bloomfilter.o ../src/markdups.o

