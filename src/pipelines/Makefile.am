EXTRA_DIST = \
	quality_assessment.sh.in \
	trimming.sh.in

r_scriptsdir = $(pkgdatadir)/r_scripts
pipelinesdir = $(pkgdatadir)/pipelines

pipelines_DATA = \
	quality_assessment.sh \
	trimming.sh

quality_assessment.sh: quality_assessment.sh.in
	sed \
		-e 's|@r_scriptsdir[@]|$(r_scriptsdir)|g' \
		-e 's|@bindir[@]|$(bindir)|g' \
		< "$<" > "$@"

trimming.sh: trimming.sh.in
	sed \
		-e 's|@bindir[@]|$(bindir)|g' \
		< "$<" > "$@"

CLEANFILES = \
	quality_assessment.sh \
	trimming.sh

MAINTAINERCLEANFILES = \
	Makefile.in
