/**
 *
 */

#ifndef __NGS_FASTA_FLEX_H__
#define __NGS_FASTA_FLEX_H__

#include "ngs_fasta.h"

/**
 * A simple and robust flex-based parser.
 */

typedef struct _FastaIterFlex FastaIterFlex;

void            iter_fasta_flex          (const char      *path,
                                          FastaIterFunc    func,
                                          void            *func_data,
                                          GError         **error);

FastaIterFlex*  fasta_iter_new_flex      (const char      *path,
                                          GError         **error);

FastaSeq*       fasta_iter_next_flex     (FastaIterFlex   *iter);
                                        
void            fasta_iter_free_flex     (FastaIterFlex   *iter);

#endif /* __NGS_FASTA_FLEX_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
