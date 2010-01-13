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

extern int   fastq_qual_char_2_int[];
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

/***********/
/* FastqDB */
/***********/

typedef struct _FastqDB FastqDB;

struct _FastqDB
{
  GHashTable   *index;

  gchar        *seqs;
  gchar        *quals;

  unsigned int  seq_size;
  unsigned int  n_seqs;

  unsigned long alloc_size;
};

FastqDB* fastq_db_new  (void);

void     fastq_db_load (FastqDB     *db,
                        const char  *path,
                        GError     **error);

void     fastq_db_free (FastqDB    *db);

/**
 * Get the option group for the fastq parsing system
 */

GOptionGroup* get_fastq_option_group (void);


#endif /* __NGS_FASTQ_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
