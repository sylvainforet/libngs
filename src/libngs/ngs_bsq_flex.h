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

/**
 * A faster parser that uses static arrays.
 * The size of the character fields is limited by the size of the static arrays.
 */
void iter_bsq_flex_ugly (char         *path,
                         BsqIterFunc   func,
                         void         *func_data,
                         GError      **error);

/**
 * A faster parser that does not do any memcopy or allocations.
 * The BsqRecord is field directly by pointer to the line.
 */
void iter_bsq_flex_line (char         *path,
                         BsqIterFunc   func,
                         void         *func_data,
                         GError      **error);

#endif /* __NGS_BSQ_FLEX_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
