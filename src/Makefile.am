EXTRA_DIST = \
	make_qual_table.py

bin_PROGRAMS = \
	fastq2fasta \
	letter_qual

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

letter_qual_SOURCES = \
	fastq.h \
	fastq.c \
	fastq_utils.h \
	fastq_utils.c \
	letter_qual.c

MAINTAINERCLEANFILES = \
	Makefile.in
