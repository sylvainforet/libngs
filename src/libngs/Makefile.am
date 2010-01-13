EXTRA_DIST = \
	make_qual_table.py

AM_CPPFLAGS = \
	@GLIB_CFLAGS@

AM_LDFLAGS = \
	@GLIB_LIBS@

lex.FlexFasta_.c: ngs_fasta_flex.l ngs_fasta_flex.h
	@LEX@ ngs_fasta_flex.l

lex.FlexFastq_.c: ngs_fastq_flex.l ngs_fastq_flex.h
	@LEX@ ngs_fastq_flex.l

lex.FlexFastqUgly_.c: ngs_fastq_flex_ugly.l ngs_fastq_flex.h
	@LEX@ ngs_fastq_flex_ugly.l

lex.FlexBsq_.c: ngs_bsq_flex.l ngs_bsq_flex.h
	@LEX@ ngs_bsq_flex.l

lex.FlexBsqLine_.c: ngs_bsq_flex_line.l ngs_bsq_flex.h
	@LEX@ ngs_bsq_flex_line.l

lex.FlexBsqUgly_.c: ngs_bsq_flex_ugly.l ngs_bsq_flex.h
	@LEX@ ngs_bsq_flex_ugly.l

noinst_LTLIBRARIES = \
	libngs.la

libngs_la_SOURCES = \
	ngs_utils.h \
	ngs_utils.c \
	ngs_fasta.h \
	ngs_fasta.c \
	ngs_fasta_flex.h \
	lex.FlexFasta_.c \
	ngs_fastq.h \
	ngs_fastq.c \
	ngs_fastq_qual_tables.c \
	ngs_fastq_flex.h \
	lex.FlexFastq_.c \
	lex.FlexFastqUgly_.c \
	ngs_bsq.h \
	ngs_bsq.c \
	ngs_bsq_flex.h \
	lex.FlexBsq_.c \
	lex.FlexBsqUgly_.c \
	lex.FlexBsqLine_.c \
	ngs_seq_db.h \
	ngs_seq_db.c

CLEANFILES = \
	lex.FlexFasta_.c \
	lex.FlexFastq_.c \
	lex.FlexFastqUgly_.c \
	lex.FlexBsq_.c \
	lex.FlexBsqUgly_.c \
	lex.FlexBsqLine_.c

MAINTAINERCLEANFILES = \
	Makefile.in
