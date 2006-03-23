AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = mincbet
mincbet_SOURCES = bet.c
mincbet_LDADD = libss.a

# lib_LIBRARIES = libss.a
noinst_LIBRARIES = libss.a

libss_a_SOURCES = \
   libss/libavw.c \
   libss/libminc.c \
   libss/libpicio.c \
   libss/libss.c \
   libss/libtessa.c

noinst_HEADERS = \
   libss/libavw.h \
   libss/libminc.h \
   libss/libpicio.h \
   libss/libss.h \
   libss/libssbase.h \
   libss/libtessa.h

m4_files = m4/mni_REQUIRE_LIB.m4 \
           m4/mni_REQUIRE_MNILIBS.m4 \
           m4/mni_REQUIRE_OPENINVENTOR.m4 \
           m4/mni_cxx_have_koenig_lookup.m4 \
           m4/smr_CGAL_MAKEFILE.m4 \
           m4/smr_OPTIONAL_LIB.m4 \
           m4/smr_REQUIRED_LIB.m4 \
           m4/smr_WITH_BUILD_PATH.m4

EXTRA_DIST = \
   $(m4_files) \
   betfunc \
   bet.tcl \
   bet_proc.tcl \
   doc/index.html \
   epm-header.in

