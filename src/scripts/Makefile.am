bin_SCRIPTS = \
	fastq_sample_mate_pairs.sh

EXTRA_DIST = \
	fastq_sample_mate_pairs.sh.in

fastq_sample_mate_pairs.sh: fastq_sample_mate_pairs.sh.in
	sed \
		-e 's|@bindir[@]|$(bindir)|g' \
		< "$<" > "$@"

CLEANFILES = \
	$(bin_SCRIPTS)

MAINTAINERCLEANFILES = \
	Makefile.in
