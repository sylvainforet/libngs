/**
 *
 */

#ifndef __FASTA_FLEX_H__
#define __FASTA_FLEX_H__

#include <glib.h>

#include "fasta.h"

/**
 * A simple and robust flex-based parser.
 */
void iter_fasta_flex (char         *path,
                      FastaIterFunc func,
                      void         *func_data,
                      GError      **error);

#endif /* __FASTA_FLEX_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
