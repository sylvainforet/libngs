/**
 *
 */

#include <ngs_utils.h>

char rev_table[128] =
{
  ['A'] = 'T',
  ['T'] = 'A',
  ['G'] = 'C',
  ['C'] = 'G',
  ['N'] = 'N'
};


char*
rev_comp_in_place (char         *seq,
                   unsigned long size)
{
  const unsigned long middle = (size + 1) / 2;
  unsigned long       i;

  for (i = 0; i < middle; i++)
    {
      const char c      = seq[i];
      seq[i]            = rev_table[(int)seq[size - i - 1]];
      seq[size - i - 1] = rev_table[(int)c];
    }
  return seq;
}

char*
rev_in_place (char         *seq,
              unsigned long size)
{
  const unsigned long middle = (size + 1) / 2;
  unsigned long       i;

  for (i = 0; i < middle; i++)
    {
      const char c      = seq[i];
      seq[i]            = seq[size - i - 1];
      seq[size - i - 1] = c;
    }
  return seq;
}

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
