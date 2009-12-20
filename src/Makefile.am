EXTRA_DIST = \
	make_qual_table.py

AM_CPPFLAGS = \
	@GLIB_CFLAGS@

AM_LDFLAGS = \
	@GLIB_LIBS@

lex.FlexFastq_.c: fastq_flex.l fastq_flex.h
	@LEX@ fastq_flex.l

lex.FlexFastqUgly_.c: fastq_flex_ugly.l fastq_flex.h
	@LEX@ fastq_flex_ugly.l

noinst_LTLIBRARIES = \
	libfastq.la

libfastq_la_SOURCES = \
	fastq.h \
	fastq.c \
	fastq_qual_tables.c \
	fastq_utils.h \
	fastq_utils.c \
	fastq_flex.h \
	lex.FlexFastq_.c \
	lex.FlexFastqUgly_.c

bin_PROGRAMS = \
	fastq2fasta \
	letter_qual

fastq2fasta_LDADD = \
	libfastq.la
fastq2fasta_SOURCES = \
	fastq2fasta.c

letter_qual_LDADD = \
	libfastq.la
letter_qual_SOURCES = \
	letter_qual.c

CLEANFILES = \
	lex.FlexFastq_.c \
	lex.FlexFastqUgly_.c

MAINTAINERCLEANFILES = \
	Makefile.in
