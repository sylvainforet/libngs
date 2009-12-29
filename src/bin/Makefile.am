AM_CPPFLAGS = \
	@GLIB_CFLAGS@

AM_LDFLAGS = \
	@GLIB_LIBS@

INCLUDES = \
	-I$(top_srcdir)/src/libngs

bin_PROGRAMS = \
	fastq2fasta \
	letter_qual \
	bsq_summary

fastq2fasta_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
fastq2fasta_SOURCES = \
	fastq2fasta.c

letter_qual_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
letter_qual_SOURCES = \
	letter_qual.c

bsq_summary_LDADD = \
	$(top_builddir)/src/libngs/libngsbsq.la
bsq_summary_SOURCES = \
	bsq_summary.c

MAINTAINERCLEANFILES = \
	Makefile.in
