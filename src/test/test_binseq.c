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

#include <stdlib.h>
#include <string.h>

#include "ngs_binseq.h"

int
main (int    argc,
      char **argv)
{
  char          *name = "seq1";
  char          *seq  = "AATT" "GCAT" "GCCA" "TTAG";
  BinSeq        *bin  = bin_seq_new (name, seq, strlen (seq));;
  char          *tmp;
  unsigned char *tmpp;

  g_print ("%s %s %ld\n", name, seq, strlen (seq));
  tmp = bin_to_char (bin->seq, bin->size);
  g_print ("%s %s %ld\n", name, tmp, strlen (tmp));
  g_free (tmp);

  tmpp = bin_revcomp_mult4 (bin->seq, bin->size);
  tmp  = bin_to_char (tmpp, bin->size);
  g_print ("%s %s %ld\n", name, tmp, strlen (tmp));
  g_free (tmp);
  g_free (tmpp);

  bin_seq_free (bin);

  return 0;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
