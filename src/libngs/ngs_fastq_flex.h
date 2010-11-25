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
void iter_fastq_flex (const char   *path,
                      FastqIterFunc func,
                      void         *func_data,
                      GError      **error);

#endif /* __NGS_FASTQ_FLEX_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
