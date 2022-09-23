AM_CXXFLAGS = \
	      -I ../libs/xxhash_cpp/include \
	      -std c++17 \
	      $(BOOST_CXXFLAGS)

AM_LDFLAGS = \
	     $(BOOST_LDFLAGS)

bin_PROGRAMS = streammd

streammd_SOURCES = bloomfilter.cxx

#LDADD = ../libs/libxxhash.a # $(INTLLIBS)
