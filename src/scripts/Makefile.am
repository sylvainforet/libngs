bin_SCRIPTS = \
	fastq_sample_mate_pairs.sh

EXTRA_DIST = \
	fastq_sample_mate_pairs.sh.in

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
