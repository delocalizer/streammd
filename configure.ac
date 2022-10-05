# Initialize autoconf
AC_INIT(streammd,
	m4_esyscmd([grep 'define STREAMMD_VERSION ' src/version.h|cut -f 2 -d '"']))

# Initialize libtool
LT_INIT

# Store auxiliary build files in their own dir
AC_CONFIG_AUX_DIR([build-aux])

# Initialize automake
AM_INIT_AUTOMAKE([subdir-objects])

# Check for C++ compiler
AC_PROG_CXX

# Source macros
AC_CONFIG_MACRO_DIR([m4])

# What to make 
AC_CONFIG_FILES([Makefile
		 src/Makefile
		 tests/Makefile])

# Generate the output
AC_OUTPUT
