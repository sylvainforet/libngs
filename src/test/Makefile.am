AM_CPPFLAGS = \
	@GLIB_CFLAGS@

AM_LDFLAGS = \
	@GLIB_LIBS@

INCLUDES = \
	-I$(top_srcdir)/src/libngs

noinst_PROGRAMS = \
	test_fasta

test_fasta_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
test_fasta_SOURCES = \
	test_fasta.c

MAINTAINERCLEANFILES = \
	Makefile.in
