bin_SCRIPTS = \
	extract_CG_contig.py \
	fastq_sample_mate_pairs.sh \
	interleave.py

EXTRA_DIST = \
	extract_CG_contig.py \
	fastq_sample_mate_pairs.sh.in \
	interleave.py \
	plot_coverage.R

fastq_sample_mate_pairs.sh: fastq_sample_mate_pairs.sh.in
	sed \
		-e 's|@bindir[@]|$(bindir)|g' \
		< "$<" > "$@"

r_scriptsdir = $(pkgdatadir)/r_scripts

r_scripts_DATA = \
	plot_coverage.R

CLEANFILES = \
	$(bin_SCRIPTS)

MAINTAINERCLEANFILES = \
	Makefile.in
