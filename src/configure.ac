AC_INIT([Vrap], [0.9], [lance@slac.stanford.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX

AC_CONFIG_FILES([
    Makefile
	path.h
])

PKG_CHECK_MODULES([LHAPDF], [lhapdf >= 6])
AM_PROG_AR

AC_PROG_RANLIB
AC_OUTPUT
