EXTRA_DIST = \
	make_qual_table.py \
	ngs_bsq_flex.l \
	ngs_fasta_flex.l \
	ngs_fastq_flex.l

AM_CPPFLAGS = \
	@GLIB_CFLAGS@

AM_LDFLAGS = \
	@GLIB_LIBS@

lex.FlexFasta_.c: ngs_fasta_flex.l ngs_fasta_flex.h
	@LEX@ ngs_fasta_flex.l

lex.FlexFastq_.c: ngs_fastq_flex.l ngs_fastq_flex.h
	@LEX@ ngs_fastq_flex.l

lex.FlexBsq_.c: ngs_bsq_flex.l ngs_bsq_flex.h
	@LEX@ ngs_bsq_flex.l

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
	ngs_bsq.h \
	ngs_bsq.c \
	ngs_bsq_flex.h \
	lex.FlexBsq_.c \
	ngs_methylation.h \
	ngs_methylation.c \
	ngs_seq_db.h \
	ngs_seq_db.c \
	ngs_binseq.h \
	ngs_binseq.c \
	ngs_binseq_tables.c \
	ngs_memalloc.h \
	ngs_memalloc.c \
	ngs_kmerhash.h \
	ngs_kmerhash.c

CLEANFILES = \
	lex.FlexFasta_.c \
	lex.FlexFastq_.c \
	lex.FlexBsq_.c

MAINTAINERCLEANFILES = \
	Makefile.in
