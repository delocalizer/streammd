# Initialize autoconf
AC_INIT(streammd, 0.1)

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

# What to make 
AC_CONFIG_FILES([
		 Makefile
		 src/Makefile])
	# when libs has source in it and needs building;
	# currently it's just got hpp file that can be statically linked
		 #libs/Makefile])

# Generate the output
AC_OUTPUT
