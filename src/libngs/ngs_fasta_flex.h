/**
 *
 */

#ifndef __NGS_FASTA_FLEX_H__
#define __NGS_FASTA_FLEX_H__

#include "ngs_fasta.h"

/**
 * A simple and robust flex-based parser.
 */
void iter_fasta_flex (char         *path,
                      FastaIterFunc func,
                      void         *func_data,
                      GError      **error);

#endif /* __NGS_FASTA_FLEX_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
