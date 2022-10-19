AM_CXXFLAGS = \
	      -Wall \
	      -Wextra \
	      -std=c++17 \
	      -I ../external

bin_PROGRAMS = streammd

streammd_SOURCES = bloomfilter.cxx markdups.cxx streammd.cxx
streammd_LDADD = ../external/xxHash/libxxhash.so

include_HEADERS = *.h
