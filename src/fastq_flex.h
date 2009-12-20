/**
 *
 */

#ifndef __FASTQ_FLEX_H__
#define __FASTQ_FLEX_H__

#include <glib.h>

#include "fastq.h"

void iter_fastq_flex (char         *path,
                      FastqIterFunc func,
                      void         *func_data,
                      GError      **error);

#endif /* __FASTQ_FLEX_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
