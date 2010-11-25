/**
 *
 */

#ifndef __NGS_FASTQ_FLEX_H__
#define __NGS_FASTQ_FLEX_H__

#include <glib.h>

#include "ngs_fastq.h"

/**
 * A simple and robust flex-based parser.
 */

typedef struct  _FastqIterFlex FastqIterFlex;

void            iter_fastq_flex         (const char    *path,
                                         FastqIterFunc  func,
                                         void          *func_data,
                                         GError       **error);

FastqIterFlex*  fastq_iter_new_flex     (const char    *path,
                                         GError       **error);

FastqSeq*       fastq_iter_next_flex    (FastqIterFlex *iter);

void            fastq_iter_free_flex    (FastqIterFlex *iter);

#endif /* __NGS_FASTQ_FLEX_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
