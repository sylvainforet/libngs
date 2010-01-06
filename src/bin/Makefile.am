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
	trim_fastq \
	sample_fastq \
	letter_qual \
	qual_summary \
	read_qual_summary \
	bsq_summary \
	bsq2gff

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

bsq2gff_LDADD = \
	$(top_builddir)/src/libngs/libngsbsq.la
bsq2gff_SOURCES = \
	bsq2gff.c

qual_summary_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
qual_summary_SOURCES = \
	qual_summary.c

read_qual_summary_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
read_qual_summary_SOURCES = \
	read_qual_summary.c

trim_fastq_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
trim_fastq_SOURCES = \
	trim_fastq.c

sample_fastq_LDADD = \
	$(top_builddir)/src/libngs/libngsfastq.la
sample_fastq_SOURCES = \
	sample_fastq.c

MAINTAINERCLEANFILES = \
	Makefile.in
