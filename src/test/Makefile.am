AM_CPPFLAGS = \
	@GLIB_CFLAGS@

AM_LDFLAGS = \
	@GLIB_LIBS@

INCLUDES = \
	-I$(top_srcdir)/src/libngs

noinst_PROGRAMS = \
	test_fasta \
	test_fasta_iter \
	test_fastq_iter \
	test_cg \
	test_binseq

test_fasta_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
test_fasta_SOURCES = \
	test_fasta.c

test_fasta_iter_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
test_fasta_iter_SOURCES = \
	test_fasta_iter.c

test_fastq_iter_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
test_fastq_iter_SOURCES = \
	test_fastq_iter.c

test_cg_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
test_cg_SOURCES = \
	test_cg.c

test_binseq_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
test_binseq_SOURCES = \
	test_binseq.c

MAINTAINERCLEANFILES = \
	Makefile.in
