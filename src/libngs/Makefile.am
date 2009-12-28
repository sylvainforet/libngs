EXTRA_DIST = \
	make_qual_table.py

AM_CPPFLAGS = \
	@GLIB_CFLAGS@

AM_LDFLAGS = \
	@GLIB_LIBS@

lex.FlexFasta_.c: fasta_flex.l fasta_flex.h
	@LEX@ fasta_flex.l

lex.FlexFastq_.c: fastq_flex.l fastq_flex.h
	@LEX@ fastq_flex.l

lex.FlexFastqUgly_.c: fastq_flex_ugly.l fastq_flex.h
	@LEX@ fastq_flex_ugly.l

lex.FlexBsq_.c: bsq_flex.l bsq_flex.h
	@LEX@ bsq_flex.l

lex.FlexBsqLine_.c: bsq_flex_line.l bsq_flex.h
	@LEX@ bsq_flex_line.l

lex.FlexBsqUgly_.c: bsq_flex_ugly.l bsq_flex.h
	@LEX@ bsq_flex_ugly.l

noinst_LTLIBRARIES = \
	libfasta.la \
	libfastq.la \
	libbsq.la

libfasta_la_SOURCES = \
	fasta.h \
	fasta.c \
	fasta_flex.h \
	lex.FlexFasta_.c

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
	lex.FlexBsqUgly_.c \
	lex.FlexBsqLine_.c

bin_PROGRAMS = \
	test_fasta \
	fastq2fasta \
	letter_qual \
	bsq_summary

test_fasta_LDADD = \
	libfasta.la
test_fasta_SOURCES = \
	test_fasta.c

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
	lex.FlexFasta_.c \
	lex.FlexFastq_.c \
	lex.FlexFastqUgly_.c \
	lex.FlexBsq_.c \
	lex.FlexBsqUgly_.c \
	lex.FlexBsqLine_.c

MAINTAINERCLEANFILES = \
	Makefile.in
