/**
 *
 */

#ifndef __NGS_FASTQ_H__
#define __NGS_FASTQ_H__

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

#define FASTQ_QUAL_0        64 /* '@' */
#define FASTQ_QUAL_0_SANGER 37 /* '%' */

extern char *fastq_qual_char_2_string[];

FastqSeq*   fastq_seq_new    (void);

void        fastq_seq_free   (FastqSeq *fastq);

FastqSeq*   fastq_seq_copy   (FastqSeq *fastq);

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

void iter_fastq (const char   *path,
                 FastqIterFunc func,
                 void         *data,
                 GError      **error);

/**
 * Get the option group for the fastq parsing system
 */

GOptionGroup* get_fastq_option_group (void);


/*************/
/* FastqIter */
/*************/

typedef struct _FastqIter FastqIter;

struct _FastqIter
{
  void *private;
};

FastqIter* fastq_iter_new  (const char *path,
                            GError    **error);

FastqSeq*  fastq_iter_next (FastqIter  *iter);

void       fastq_iter_free (FastqIter  *iter);


#endif /* __NGS_FASTQ_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
