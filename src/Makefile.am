AM_CXXFLAGS = \
	      -Wall \
	      -Wextra \
	      -std=c++17 \
	      -I ../external 

bin_PROGRAMS = streammd

streammd_SOURCES = streammd.cxx markdups.cxx bloomfilter.cxx
#streammd_LDADD = libbloomfilter.a libmarkdups.a

#noinst_LIBRARIES = libbloomfilter.a libmarkdups.a
#libbloomfilter_a_SOURCES = bloomfilter.cxx
#libmarkdups_a_SOURCES = markdups.cxx

include_HEADERS = *.h
