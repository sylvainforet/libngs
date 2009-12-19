/**
 *
 */

#ifndef __FASTQ_H__
#define __FASTQ_H__

#include <glib.h>

/************/
/* FastqSeq */
/************/

typedef struct _FastqSeq FastqSeq;

struct _FastqSeq
{
  char *name;
  char *seq;
  char *qual;
  int   size;
};

extern int   fastq_qual_char_2_int[];
extern char *fastq_qual_char_2_string[];

FastqSeq*   fastq_seq_new    (void);

void        fastq_seq_free   (FastqSeq *fastq);

/**
 * Signature for the functions to be called by iter_fastq.
 * If the function returns 0, the iteration is interupted,
 * otherwise, the iteration continues.
 */

typedef int (*FastqIterFunc) (FastqSeq *fastq,
                              void     *data);

/**
 * Iterates over all FastqSeq in a file.
 * If path is '-', reads from stdin.
 */

void iter_fastq (char         *path,
                 FastqIterFunc func,
                 void         *data,
                 GError      **error);


#endif /* __FASTQ_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
