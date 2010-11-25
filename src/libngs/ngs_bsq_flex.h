/**
 *
 */

#ifndef __NGS_BSQ_FLEX_H__
#define __NGS_BSQ_FLEX_H__

#include <glib.h>

#include "ngs_bsq.h"

/**
 * A simple and robust flex-based parser.
 */
void iter_bsq_flex (char         *path,
                    BsqIterFunc   func,
                    void         *func_data,
                    GError      **error);

#endif /* __NGS_BSQ_FLEX_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
