#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([smo_svmpack], [alpha], [patjnichols@gmail.com])
AC_CONFIG_SRCDIR([svm_options.cpp])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG([C++])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h omp.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_CHECK_FUNCS([clock_gettime memset memcpy rint setenv sqrt posix_memalign malloc align_malloc])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
