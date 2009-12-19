EXTRA_DIST = \
	make_qual_table.py

bin_PROGRAMS = \
	fastq2fasta

AM_CPPFLAGS = \
	@GLIB_CFLAGS@

AM_LDFLAGS = \
	@GLIB_LIBS@

fastq2fasta_SOURCES = \
	fastq.h \
	fastq.c \
	fastq_utils.h \
	fastq_utils.c \
	fastq2fasta.c

MAINTAINERCLEANFILES = \
	Makefile.in
