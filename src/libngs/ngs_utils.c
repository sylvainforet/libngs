/* Copyright (C) 2010  Sylvain FORET
 *
 * This file is part of libngs.
 *
 * libngs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *                                                                       
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *                                                                       
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
