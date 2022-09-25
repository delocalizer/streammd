AM_CXXFLAGS = \
	      -std=c++17 \
	      -I ../external \
	      $(BOOST_CXXFLAGS)

AM_LDFLAGS = \
	     $(BOOST_LDFLAGS)

bin_PROGRAMS = streammd

streammd_SOURCES = markdups.cxx
streammd_LDADD = libbloomfilter.a

noinst_LIBRARIES = libbloomfilter.a
libbloomfilter_a_SOURCES = bloomfilter.cxx

