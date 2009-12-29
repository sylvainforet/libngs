/**
 *
 */

#include <ngs_utils.h>

GQuark
ngs_error_quark (void)
{
  static GQuark error_quark = 0;
  if (error_quark == 0)
    return g_quark_from_static_string ("ngs-error-quark");
  return error_quark;
}


/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
