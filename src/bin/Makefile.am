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
	fastq_interleave \
	fastq_letter_qual \
	fastq_letter_pos \
	fastq_qual_length_summary \
	fastq_read_qual_summary \
	fastq_revcomp \
	fastq_sample \
	fastq_split \
	fastq_trim \
	bsq2gff \
	bsq_coverage_distribution \
	bsq_empty_regions_content \
	bsq_methylation_counts \
	bsq_summary \
	cg_fetch \
	cg_merge \
	cg_meth_count \
	cg_meth_dist \
	kmers_count \
	kmers_count_array16 \
	kmers_bin_tool

fastq2fasta_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq2fasta_SOURCES = \
	fastq2fasta.c

fastq_base_qual_summary_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_base_qual_summary_SOURCES = \
	fastq_base_qual_summary.c

fastq_interleave_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_interleave_SOURCES = \
	fastq_interleave.c

fastq_fetch_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_fetch_SOURCES = \
	fastq_fetch.c

fastq_letter_pos_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_letter_pos_SOURCES = \
	fastq_letter_pos.c

fastq_letter_qual_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_letter_qual_SOURCES = \
	fastq_letter_qual.c

fastq_qual_length_summary_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_qual_length_summary_SOURCES = \
	fastq_qual_length_summary.c

fastq_read_qual_summary_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_read_qual_summary_SOURCES = \
	fastq_read_qual_summary.c

fastq_revcomp_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_revcomp_SOURCES = \
	fastq_revcomp.c

fastq_sample_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_sample_SOURCES = \
	fastq_sample.c

fastq_split_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
fastq_split_SOURCES = \
	fastq_split.c

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

bsq_methylation_counts_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
bsq_methylation_counts_SOURCES = \
	bsq_methylation_counts.c

bsq_summary_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
bsq_summary_SOURCES = \
	bsq_summary.c

cg_fetch_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
cg_fetch_SOURCES = \
	cg_fetch.c

cg_merge_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
cg_merge_SOURCES = \
	cg_merge.c

cg_meth_count_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
cg_meth_count_SOURCES = \
	cg_meth_count.c

cg_meth_dist_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
cg_meth_dist_SOURCES = \
	cg_meth_dist.c

kmers_count_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
kmers_count_SOURCES = \
	kmers_count.c

kmers_count_array16_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
kmers_count_array16_SOURCES = \
	kmers_count_array16.c

kmers_bin_tool_LDADD = \
	$(top_builddir)/src/libngs/libngs.la
kmers_bin_tool_SOURCES = \
	kmers_bin_tool.c

MAINTAINERCLEANFILES = \
	Makefile.in
