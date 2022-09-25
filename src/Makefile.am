AM_CXXFLAGS = \
	      $(BOOST_CXXFLAGS) \
	      -DBOOST_LOG_DYN_LINK \
	      -Wall \
	      -Wextra \
	      -std=c++17 \
	      -I ../external 

AM_LDFLAGS = \
	     -lboost_log \
	     -lboost_log_setup \
	     $(BOOST_LDFLAGS)

bin_PROGRAMS = streammd

streammd_SOURCES = markdups.cxx
streammd_LDADD = libbloomfilter.a

noinst_LIBRARIES = libbloomfilter.a
libbloomfilter_a_SOURCES = bloomfilter.cxx

