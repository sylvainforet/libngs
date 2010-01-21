bin_SCRIPTS = \
	extract_CG_contig.py \
	fastq_sample_mate_pairs.sh \
	interleave.py \
	table2fastq.py

EXTRA_DIST = \
	extract_CG_contig.py \
	fastq_sample_mate_pairs.sh.in \
	interleave.py \
	plot_coverage.R \
	table2fastq.py

fastq_sample_mate_pairs.sh: fastq_sample_mate_pairs.sh.in
	sed \
		-e 's|@bindir[@]|$(bindir)|g' \
		< "$<" > "$@"

r_scriptsdir = $(pkgdatadir)/r_scripts

r_scripts_DATA = \
	plot_coverage.R

CLEANFILES = \
	fastq_sample_mate_pairs.sh

MAINTAINERCLEANFILES = \
	Makefile.in
