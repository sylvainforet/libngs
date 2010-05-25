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
  char   *name = "seq1";
  char   *seq  = "AATT" "GCAT" "GCCA" "TT";
  BinSeq *bin  = bin_seq_new (name, seq, strlen (seq));;
  char   *tmp;

  g_print ("%s %s %ld\n", name, seq, strlen (seq));
  tmp = bin_to_char (bin->seq, bin->size);
  g_print ("%s %s %ld\n", name, tmp, strlen (tmp));

  g_free (tmp);
  bin_seq_free (bin);

  return 0;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
