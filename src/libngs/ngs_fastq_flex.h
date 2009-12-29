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
void iter_fastq_flex (char         *path,
                      FastqIterFunc func,
                      void         *func_data,
                      GError      **error);

/**
 * A faster parser that uses static arrays.
 * The size of the FastqSeq fields is limited by the size of the static arrays.
 */
void iter_fastq_flex_ugly (char         *path,
                           FastqIterFunc func,
                           void         *func_data,
                           GError      **error);

#endif /* __NGS_FASTQ_FLEX_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
