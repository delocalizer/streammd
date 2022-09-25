# Initialize autoconf
AC_INIT(streammd,
	m4_esyscmd([grep 'define STREAMMD_VERSION ' src/version.h|cut -f 2 -d '"']))

# Initialize libtool
LT_INIT

# Store auxiliary build files in their own dir
AC_CONFIG_AUX_DIR([build-aux])

# Initialize automake
AM_INIT_AUTOMAKE

# Check for C++ compiler
AC_PROG_CXX

# Source macros
AC_CONFIG_MACRO_DIR([m4])

# Check for Boost
AX_BOOST_BASE([1.56],
	      [],
	      [AC_MSG_ERROR([streammd needs Boost >= 1.56 but it was not found in your system])])
AX_BOOST_SYSTEM
AX_BOOST_LOG
AX_BOOST_LOG_SETUP

# What to make 
AC_CONFIG_FILES([Makefile
		 src/Makefile])

# Generate the output
AC_OUTPUT
