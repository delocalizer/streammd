AM_CXXFLAGS = \
	      $(BOOST_CXXFLAGS) \
	      -Wall \
	      -Wextra \
	      -std=c++17 \
	      -I ../external 

AM_LDFLAGS = \
	     $(BOOST_LDFLAGS)

bin_PROGRAMS = streammd

streammd_SOURCES = markdups.cxx
streammd_LDADD = libbloomfilter.a

noinst_LIBRARIES = libbloomfilter.a
libbloomfilter_a_SOURCES = bloomfilter.cxx

