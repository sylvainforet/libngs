EXTRA_DIST = \
	make_qual_table.py

bin_PROGRAMS = \
	fastq2fasta \
	letter_qual

AM_CPPFLAGS = \
	@GLIB_CFLAGS@

AM_LDFLAGS = \
	@GLIB_LIBS@

lex.FlexFastq_.c: fastq_flex.l fastq_flex.h
	@LEX@ fastq_flex.l

fastq2fasta_SOURCES = \
	fastq.h \
	fastq.c \
	lex.FlexFastq_.c \
	fastq_flex.h \
	fastq_utils.h \
	fastq_utils.c \
	fastq2fasta.c

letter_qual_SOURCES = \
	fastq.h \
	fastq.c \
	lex.FlexFastq_.c \
	fastq_flex.h \
	fastq_utils.h \
	fastq_utils.c \
	letter_qual.c

CLEANFILES = \
	lex.FlexFastq_.c

MAINTAINERCLEANFILES = \
	Makefile.in
