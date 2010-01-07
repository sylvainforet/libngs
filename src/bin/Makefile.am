EXTRA_DIST = \
	template.c

AM_CPPFLAGS = \
	@GLIB_CFLAGS@

AM_LDFLAGS = \
	@GLIB_LIBS@

INCLUDES = \
	-I$(top_srcdir)/src/libngs

bin_PROGRAMS = \
	fastq2fasta \
	fastq_base_qual_summary \
	fastq_fetch \
	fastq_letter_qual \
	fastq_read_qual_summary \
	fastq_sample \
	fastq_trim \
	bsq_summary \
	bsq2gff

fastq2fasta_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
fastq2fasta_SOURCES = \
	fastq2fasta.c

fastq_letter_qual_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
fastq_letter_qual_SOURCES = \
	fastq_letter_qual.c

bsq_summary_LDADD = \
	$(top_builddir)/src/libngs/libngsbsq.la
bsq_summary_SOURCES = \
	bsq_summary.c

bsq2gff_LDADD = \
	$(top_builddir)/src/libngs/libngsbsq.la
bsq2gff_SOURCES = \
	bsq2gff.c

fastq_base_qual_summary_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
fastq_base_qual_summary_SOURCES = \
	fastq_base_qual_summary.c

fastq_read_qual_summary_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
fastq_read_qual_summary_SOURCES = \
	fastq_read_qual_summary.c

fastq_trim_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
fastq_trim_SOURCES = \
	fastq_trim.c

fastq_sample_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
fastq_sample_SOURCES = \
	fastq_sample.c

fastq_fetch_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
fastq_fetch_SOURCES = \
	fastq_fetch.c

MAINTAINERCLEANFILES = \
	Makefile.in
