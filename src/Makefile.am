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

lex.FlexBsq_.c: bsq_flex.l bsq_flex.h
	@LEX@ bsq_flex.l

lex.FlexBsqUgly_.c: bsq_flex_ugly.l bsq_flex.h
	@LEX@ bsq_flex_ugly.l

noinst_LTLIBRARIES = \
	libfastq.la \
	libbsq.la

libfastq_la_SOURCES = \
	fastq.h \
	fastq.c \
	fastq_qual_tables.c \
	fastq_utils.h \
	fastq_utils.c \
	fastq_flex.h \
	lex.FlexFastq_.c \
	lex.FlexFastqUgly_.c

libbsq_la_SOURCES = \
	bsq.h \
	bsq.c \
	bsq_flex.h \
	lex.FlexBsq_.c \
	lex.FlexBsqUgly_.c

bin_PROGRAMS = \
	fastq2fasta \
	letter_qual \
	bsq_summary

fastq2fasta_LDADD = \
	libfastq.la
fastq2fasta_SOURCES = \
	fastq2fasta.c

letter_qual_LDADD = \
	libfastq.la
letter_qual_SOURCES = \
	letter_qual.c

bsq_summary_LDADD = \
	libbsq.la
bsq_summary_SOURCES = \
	bsq_summary.c

CLEANFILES = \
	lex.FlexFastq_.c \
	lex.FlexFastqUgly_.c \
	lex.FlexBsq_.c
	lex.FlexBsqUgly_.c

MAINTAINERCLEANFILES = \
	Makefile.in
