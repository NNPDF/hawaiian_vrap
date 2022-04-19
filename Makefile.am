
AM_CPPFLAGS = $(shell lhapdf-config --cppflags) 
AM_LDFLAGS = $(shell lhapdf-config --ldflags) 

lib_LIBRARIES = libVrap.a 
bin_PROGRAMS = Vrap Vtot 

libVrap_a_SOURCES = \
	a02.C  a02.h   \
	vegas.C     \
	dilog.C     \
	totfns.C    \
	pdf.C pdf.h      \
	CClasses.h CClasses.C    \
	DClasses.h DClasses.C    \
	IClasses.h IClasses.C    \
	NAClasses.h  vegas.h dilog.h \
	totfns.h \
	complex.h complex.C    \
	thbasics.h thbasics.C    \
	random.h random.C    \
	QCDbasics.h QCDbasics.C    \
	hardfns.h hardfns.C     \
	Vlumifns.h Vlumifns.C    \
	BornNLOfns.h BornNLOfns.C    \
	qbarqfns.h qbarqfns.C    \
	qgfns.h qgfns.C     \
	ggfns.h ggfns.C     \
	qiqjfns.h qiqjfns.C    \
	Vlumifns_LHApdf.C Vlumifns_LHApdf.h \
	LHApdf.cpp LHApdf.h EWparams.h \
	mode.C mode.h integration.h integration_tot.h \
	options.h options.hpp options.C

Vrap_SOURCES = Vrap.C
Vtot_SOURCES = Vtot.C

Vrap_LDADD = libVrap.a 
Vtot_LDADD = libVrap.a 

dist_data = INSTALL