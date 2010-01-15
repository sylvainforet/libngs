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
	bsq2gff \
	bsq_coverage_distribution \
	bsq_empty_regions_content \
	bsq_methylation_ratio \
	bsq_summary \
	cg_merge

fastq2fasta_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq2fasta_SOURCES = \
	fastq2fasta.c

fastq_base_qual_summary_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_base_qual_summary_SOURCES = \
	fastq_base_qual_summary.c

fastq_fetch_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_fetch_SOURCES = \
	fastq_fetch.c

fastq_letter_qual_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_letter_qual_SOURCES = \
	fastq_letter_qual.c

fastq_read_qual_summary_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_read_qual_summary_SOURCES = \
	fastq_read_qual_summary.c

fastq_sample_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_sample_SOURCES = \
	fastq_sample.c

fastq_trim_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_trim_SOURCES = \
	fastq_trim.c

bsq2gff_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
bsq2gff_SOURCES = \
	bsq2gff.c

bsq_coverage_distribution_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
bsq_coverage_distribution_SOURCES = \
	bsq_coverage_distribution.c

bsq_empty_regions_content_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
bsq_empty_regions_content_SOURCES = \
	bsq_empty_regions_content.c

bsq_methylation_ratio_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
bsq_methylation_ratio_SOURCES = \
	bsq_methylation_ratio.c

bsq_summary_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
bsq_summary_SOURCES = \
	bsq_summary.c

cg_merge_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
cg_merge_SOURCES = \
	cg_merge.c

MAINTAINERCLEANFILES = \
	Makefile.in
