/**
 *
 */

#ifndef __NGS_FASTA_H__
#define __NGS_FASTA_H__

#include <glib.h>

/************/
/* FastaSeq */
/************/

typedef struct _FastaSeq FastaSeq;

struct _FastaSeq
{
  char              *name;
  char              *seq;
  unsigned long int  size;
};

FastaSeq*   fasta_seq_new    (void);

void        fasta_seq_free   (FastaSeq *fasta);

/**
 * Signature for the functions to be called by iter_fasta.
 * If the function returns 0, the iteration is interupted,
 * otherwise, the iteration continues.
 */

typedef int (*FastaIterFunc) (FastaSeq *fasta,
                              void     *data);

/**
 * Iterates over all FastaSeq in a file.
 * If path is '-', reads from stdin.
 */

void iter_fasta (char         *path,
                 FastaIterFunc func,
                 void         *data,
                 GError      **error);

/**
 * Loads a fasta file into a hash table
 */

GHashTable* load_fasta_dict (char *path);

/**
 * Get the option group for the fasta parsing system
 */

GOptionGroup* get_fasta_option_group (void);

#endif /* __NGS_FASTA_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
