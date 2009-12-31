/**
 *
 */

#ifndef __FASTQ_UTILS_H__
#define __FASTQ_UTILS_H__

#include <glib.h>

typedef enum
{
  NGS_UNKNOWN_ERROR,
  NGS_PARSE_ERROR,
  NGS_IO_ERROR
}
NgsError;

GQuark ngs_error_quark (void);

#define NGS_ERROR ngs_error_quark ()

#endif /* __FASTQ_UTILS_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */