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


typedef enum
{
  NUC_A  = 0,
  NUC_C  = 1,
  NUC_G  = 2,
  NUC_T  = 3,
  NUC_NB = 4
}
Nucleotide;

unsigned char char_to_bin_table[128] =
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

char bin_to_char_table[NUC_NB] =
{
  [NUC_A] = 'A',
  [NUC_C] = 'C',
  [NUC_G] = 'G',
  [NUC_T] = 'T'
};

BinSeq*
bin_seq_new  (const char        *name,
              const char        *seq,
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

unsigned char*
char_to_bin (const char        *seq,
             unsigned long int  size)
{
  size_t         byte_size;
  unsigned char *new_seq;

  byte_size = (size + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;
  new_seq   = g_malloc0 (byte_size);
  new_seq   = char_to_bin_prealloc (new_seq, seq, size);

  return new_seq;
}

unsigned char*
char_to_bin_prealloc (unsigned char     *dest,
                      const char        *src,
                      unsigned long int  size)
{
  const unsigned long int mod_size = (size / NUCS_PER_BYTE) * NUCS_PER_BYTE;
  const unsigned char    *usrc = (unsigned char *)src;
  unsigned long int       i;
  unsigned long int       j;

  for (i = 0, j = 0; i < mod_size; i += 4, j++)
    {
      const unsigned long int i1 = i + 1;
      const unsigned long int i2 = i + 2;
      const unsigned long int i3 = i + 3;
      dest[j]  = char_to_bin_table[usrc[i]];
      dest[j] |= char_to_bin_table[usrc[i1]] << 2;
      dest[j] |= char_to_bin_table[usrc[i2]] << 4;
      dest[j] |= char_to_bin_table[usrc[i3]] << 6;
    }
  if (i < size)
    {
      unsigned int offset = BITS_PER_NUC;
      dest[j]             = char_to_bin_table[usrc[i]];
      for (++i ; i < size; i++)
        {
          dest[j] |= char_to_bin_table[usrc[i]] << offset;
          offset  += BITS_PER_NUC;
        }
    }
  return dest;
}

char*
bin_to_char (const unsigned char *seq,
             unsigned long int    size)
{
  char *new_seq;

  /* '\0' terminated */
  new_seq = g_malloc (size + 1);
  new_seq = bin_to_char_prealloc (new_seq, seq, size);
  new_seq[size] = '\0';

  return new_seq;
}

char*
bin_to_char_prealloc  (char                *dest,
                       const unsigned char *src,
                       unsigned long int    size)
{
  const unsigned long int mod_size = (size / NUCS_PER_BYTE) * NUCS_PER_BYTE;
  unsigned long int       i;
  unsigned long int       j;
  unsigned int            offset;

  for (i = 0, j = 0; i < mod_size; i += 4, j++)
    {
      const unsigned long int i1 = i + 1;
      const unsigned long int i2 = i + 2;
      const unsigned long int i3 = i + 3;
      dest[i]  = (bin_to_char_table[src[j]      & 3]);
      dest[i1] = (bin_to_char_table[src[j] >> 2 & 3]);
      dest[i2] = (bin_to_char_table[src[j] >> 4 & 3]);
      dest[i3] = (bin_to_char_table[src[j] >> 6 & 3]);
    }

  for (offset = 0; i < size; i++)
    {
      const int val = (src[j] >> offset) & 3;
      dest[i]       = bin_to_char_table[val];
      offset       += BITS_PER_NUC;
    }
  return dest;
}

unsigned char*
bin_revcomp_mult4 (const unsigned char *seq,
                   unsigned long int    size)
{
  unsigned char *revcomp;

  revcomp = g_malloc0 (size / NUCS_PER_BYTE);
  revcomp = bin_revcomp_prealloc_mult4 (revcomp, seq, size);

  return revcomp;
}

unsigned char*
bin_revcomp_prealloc_mult4 (unsigned char       *dest,
                            const unsigned char *src,
                            unsigned long int    size)
{
  unsigned long int       i;
  const unsigned long int n_bytes = size / NUCS_PER_BYTE;

  for (i = 0; i < n_bytes; i++)
    dest[n_bytes - i - 1] = revcomp_table[src[i]];

  return dest;
}

unsigned char*
bin_revcomp_inplace_mult4 (unsigned char     *seq,
                           unsigned long int  size)
{
  unsigned long int       i;
  const unsigned long int maxi = (size / NUCS_PER_BYTE + 1) / 2;

  for (i = 0; i < maxi; i++)
    {
      const unsigned long int comp_idx = size - i - 1;
      const unsigned char     tmp      = revcomp_table[seq[i]];

      seq[i]        = revcomp_table[seq[comp_idx]];
      seq[comp_idx] = tmp;
    }

  return seq;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
