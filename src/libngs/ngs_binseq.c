/**
 *
 */

#include <stdlib.h>
#include <string.h>

#include "ngs_binseq.h"


static char char_to_bin_table[128] =
{
  ['a'] = NUC_A,
  ['c'] = NUC_C,
  ['g'] = NUC_G,
  ['t'] = NUC_T,
  ['A'] = NUC_A,
  ['C'] = NUC_C,
  ['G'] = NUC_G,
  ['T'] = NUC_T,
};

static char bin_to_char_table[NUC_NB] =
{
  [NUC_A] = 'A',
  [NUC_C] = 'C',
  [NUC_G] = 'G',
  [NUC_T] = 'T'
};

BinSeq*
bin_seq_new  (char              *name,
              char              *seq,
              unsigned long int  size)
{
  BinSeq *bin;

  bin       = g_slice_new0   (BinSeq);
  bin->name = g_strdup (name);
  bin->size = size;
  bin->seq  = char_to_bin (seq, size);

  return bin;
}

void
bin_seq_free (BinSeq *bin)
{
  if (bin)
    {
      if (bin->name)
        g_free (bin->name);
      if (bin->seq)
        g_free (bin->seq);
      g_slice_free (BinSeq, bin);
    }
}

char*
char_to_bin (char              *seq,
             unsigned long int  size)
{
  char             *new_seq;
  size_t            byte_size;
  unsigned long int i;

  byte_size = (size + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;
  new_seq   = g_malloc0 (byte_size);

  for (i = 0; i < size; i++)
    {
      const unsigned long int address = i / NUCS_PER_BYTE;
      const unsigned int      offset  = (i % NUCS_PER_BYTE) * BITS_PER_NUC;
      new_seq[address]               |= char_to_bin_table[(int)seq[i]] << offset;
    }

  return new_seq;
}

char*
bin_to_char (char              *seq,
             unsigned long int  size)
{
  char             *new_seq;
  unsigned long int i;

  /* '\0' terminated */
  new_seq = g_malloc0 (size + 1);

  for (i = 0; i < size; i++)
    {
      const unsigned long int address = i / NUCS_PER_BYTE;
      const unsigned int      offset  = (i % NUCS_PER_BYTE) * BITS_PER_NUC;
      const int               val     = (seq[address] >> offset) & 3;
      new_seq[i]                      = bin_to_char_table[val];
    }

  return new_seq;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
