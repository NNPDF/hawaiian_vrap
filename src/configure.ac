AC_INIT([Vrap], [1.6], [juacrumar@lairen.eu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX

AC_CONFIG_FILES([
    Makefile
	path.h
])

PKG_CHECK_MODULES([LHAPDF], [lhapdf >= 6])
PKG_CHECK_MODULES([pineappl_capi], [pineappl_capi])
AM_PROG_AR

AC_PROG_RANLIB
AC_OUTPUT
