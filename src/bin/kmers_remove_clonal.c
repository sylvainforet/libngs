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
#include <unistd.h>

#include "ngs_binseq.h"
#include "ngs_fastq.h"
#include "ngs_kmerhash.h"
#include "ngs_utils.h"

typedef struct _FastqPair FastqPair;

struct _FastqPair
{
  FastqSeq *seq1;
  FastqSeq *seq2;
};

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char         *input_path1;
  char         *input_path2;
  char         *output_path;
  char         *output_path1;
  char         *output_path2;

  char        **coords_str;
  unsigned int *start_coords;
  unsigned int *end_coords;
  int           n_coords;

  int           verbose;
};

static void             parse_args      (CallbackData   *data,
                                         int            *argc,
                                         char         ***argv);

static void             parse_coords    (CallbackData   *data);

static int              check_sanity    (CallbackData   *data,
                                         FastqSeq       *seq1,
                                         FastqSeq       *seq2,
                                         unsigned int    start,
                                         unsigned int    end);

static KmerHashTable*   load_sequences  (CallbackData   *data);

static KmerHashTable*   remove_clonal   (CallbackData   *data,
                                         KmerHashTable  *hash_table,
                                         int             round);

static int              fastq_pair_cmp  (CallbackData   *data,
                                         FastqSeq       *seq11,
                                         FastqSeq       *seq12,
                                         FastqSeq       *seq21,
                                         FastqSeq       *seq22);

static void             write_sequences (CallbackData   *data,
                                         KmerHashTable  *hash_table);

static void             write_fastq     (FastqSeq       *fastq,
                                         GString        *buffer,
                                         GIOChannel     *channel);

int
main (int    argc,
      char **argv)
{
  CallbackData   data;
  KmerHashTable *hash_table;
  int            i;

  parse_args (&data, &argc, &argv);
  hash_table = load_sequences (&data);
  for (i = 1; i < data.n_coords; i++)
    hash_table = remove_clonal (&data, hash_table, i);

  write_sequences (&data, hash_table);
  kmer_hash_table_destroy (hash_table);

  return 0;
}

static void
parse_args (CallbackData   *data,
            int            *argc,
            char         ***argv)
{
  GOptionEntry entries[] =
    {
      {"output",  'o', 0, G_OPTION_ARG_FILENAME,     &data->output_path, "Output path prefix", NULL},
      {"coords",  'c', 0, G_OPTION_ARG_STRING_ARRAY, &data->coords_str,  "Coordinates for filtering (1-based, inclusive)", "START,END"},
      {"verbose", 'v', 0, G_OPTION_ARG_NONE,         &data->verbose,     "Verbose output", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->output_path    = NULL;
  data->output_path1   = NULL;
  data->output_path2   = NULL;
  data->coords_str     = NULL;
  data->start_coords   = NULL;
  data->end_coords     = NULL;
  data->verbose        = 0;

  context = g_option_context_new ("FILE1 FILE2 - filters out clonal sequences from a pair fastq files");
  g_option_context_add_group (context, get_fastq_option_group ());
  g_option_context_add_main_entries (context, entries, NULL);
  if (!g_option_context_parse (context, argc, argv, &error))
    {
      g_printerr ("[ERROR] Option parsing failed: %s\n", error->message);
      exit (1);
    }
  g_option_context_free (context);

  if (*argc < 3)
    {
      g_printerr ("[ERROR] No input file provided\n");
      exit (1);
    }
  data->input_path1 = (*argv)[1];
  data->input_path2 = (*argv)[2];

  parse_coords (data);

  if (data->output_path)
    {
      data->output_path1 = g_strconcat (data->output_path, ".1", NULL);
      data->output_path2 = g_strconcat (data->output_path, ".2", NULL);
    }
  else
    data->output_path = g_strdup ("-");
}

static void
write_sequences (CallbackData   *data,
                 KmerHashTable  *hash_table)
{
  KmerHashTableIter iter;
  KmerHashNode     *node;
  GIOChannel       *output_channel1;
  GIOChannel       *output_channel2;
  GString          *buffer;
  GError           *error = NULL;
  int               use_stdout = 0;

  if (data->verbose)
    g_printerr ("Writing sequences\n");

  if (data->output_path[0] == '-' && data->output_path[1] == '\0')
    {
      use_stdout     = 1;
      output_channel1 = g_io_channel_unix_new (STDOUT_FILENO);
      output_channel2 = output_channel1;
    }
  else
    {
      output_channel1 = g_io_channel_new_file (data->output_path1,
                                              "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening sequence output file failed: %s\n",
                      error->message);
          exit (1);
        }
      output_channel2 = g_io_channel_new_file (data->output_path2,
                                              "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening sequence output file failed: %s\n",
                      error->message);
          exit (1);
        }
      g_io_channel_set_encoding (output_channel2, NULL, NULL);
      g_io_channel_set_buffered (output_channel2, TRUE);
      g_io_channel_set_buffer_size (output_channel2, 4096 * 128);
    }
  g_io_channel_set_encoding (output_channel1, NULL, NULL);
  g_io_channel_set_buffered (output_channel1, TRUE);
  g_io_channel_set_buffer_size (output_channel1, 4096 * 128);

  buffer = g_string_sized_new (1024);
  kmer_hash_table_iter_init (&iter, hash_table);
  while ((node = kmer_hash_table_iter_next (&iter)) != NULL)
    {
      FastqPair *pair;

      pair = (FastqPair*) node->value;
      write_fastq (pair->seq1, buffer, output_channel1);
      write_fastq (pair->seq2, buffer, output_channel2);
      fastq_seq_free (pair->seq1);
      fastq_seq_free (pair->seq2);
      g_slice_free (FastqPair, pair);
    }

  /* Cleanup */
  g_string_free (buffer, TRUE);
  if (!use_stdout)
    {
      g_io_channel_shutdown (output_channel1, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing output file failed: %s\n", error->message);
          g_error_free (error);
          error = NULL;
        }
      g_io_channel_shutdown (output_channel2, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing output file failed: %s\n", error->message);
          g_error_free (error);
          error = NULL;
        }
      g_io_channel_unref (output_channel2);
      g_free (data->output_path1);
      g_free (data->output_path2);
    }
  g_io_channel_unref (output_channel1);
  g_free (data->output_path);
  g_free (data->start_coords);
  g_free (data->end_coords);
}

static void
write_fastq (FastqSeq   *fastq,
             GString    *buffer,
             GIOChannel *channel)
{
  GError  *error = NULL;

  g_string_truncate (buffer, 0);

  /* Sequence */
  buffer = g_string_append_c (buffer, '@');
  buffer = g_string_append (buffer, fastq->name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append (buffer, fastq->seq);
  buffer = g_string_append_c (buffer, '\n');

  /* Quality */
  buffer = g_string_append_c (buffer, '+');
  buffer = g_string_append (buffer, fastq->name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append (buffer, fastq->qual);
  buffer = g_string_append_c (buffer, '\n');

  g_io_channel_write_chars (channel,
                            buffer->str,
                            -1,
                            NULL,
                            &error);
  if (error)
    {
      g_printerr ("[ERROR] Writing sequence failed: %s\n", error->message);
      exit (1);
    }
}

static void
parse_coords (CallbackData *data)
{
  int i;

  if (data->coords_str == NULL)
    {
      g_printerr ("[ERROR] No coordinates provided, see the `-c/--coords' switch\n");
      exit (1);
    }

  for (data->n_coords = 0;
       data->coords_str[data->n_coords] != NULL;
       data->n_coords++);

  data->start_coords = g_malloc0 (data->n_coords * sizeof (*data->start_coords));
  data->end_coords   = g_malloc0 (data->n_coords * sizeof (*data->end_coords));

  for (i = 0; i < data->n_coords; i++)
    {
      char        **fields;
      unsigned int  start;
      unsigned int  end;
      int           j;

      fields = g_strsplit (data->coords_str[i], ",", 0);

      for (j = 0;
           fields[j] != NULL;
           j++);
      if (j != 2)
        {
          g_printerr ("[ERROR] Could not parse coordinates: %s\n",
                      data->coords_str[i]);
          exit (1);
        }
      start = g_ascii_strtoull (fields[0], NULL, 10);
      end   = g_ascii_strtoull (fields[1], NULL, 10);
      if (start == 0)
        {
          start = 1;
          g_printerr ("[WARNING] Start is 0, changed to 1\n");
        }
      if (end == 0)
        {
          end = 1;
          g_printerr ("[WARNING] End is 0, changed to 1\n");
        }
      if (start <= end)
        {
          data->start_coords[i] = start;
          data->end_coords[i]   = end;
        }
      else
        {
          g_printerr ("[WARNING] Start (%u) is greater than end (%d), inverting\n",
                      start, end);
          data->start_coords[i] = end;
          data->end_coords[i]   = start;
        }
      g_strfreev (fields);
    }
  g_strfreev (data->coords_str);
}

static int
check_sanity (CallbackData   *data,
              FastqSeq       *seq1,
              FastqSeq       *seq2,
              unsigned int    start,
              unsigned int    end)
{
  int i;

  /* Check size */
  if (seq1->size < end || seq2->size < end)
    return 0;

  /* Dump reads with Ns */
  for (i = start - 1; i < end; i++)
    {
      if (seq1->seq[i] == 'N' ||
          seq2->seq[i] == 'N' ||
          seq1->seq[i] == 'n' ||
          seq2->seq[i] == 'n')
        break;
    }
  if (i != end)
    return 0;

  return 1;
}

static KmerHashTable*
load_sequences  (CallbackData   *data)
{
  GError        *error = NULL;
  FastqSeq      *seq1;
  FastqSeq      *seq2;
  FastqIter     *iter1;
  FastqIter     *iter2;
  KmerHashTable *hash_table;
  unsigned char *kmer;
  unsigned int   start;
  unsigned int   end;
  unsigned int   k;
  unsigned int   k_bytes;

  if (data->verbose)
    g_printerr ("Loading sequences (round 1)\n");

  start      = data->start_coords[0];
  end        = data->end_coords[0];
  k          = end - start + 1;
  k_bytes    = (k + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;
  hash_table = kmer_hash_table_new (k_bytes * NUCS_PER_BYTE * 2);
  kmer       = g_malloc0 (2 * k_bytes * sizeof (*kmer));

  iter1 = fastq_iter_new (data->input_path1, &error);
  if (error)
    {
      g_printerr ("[ERROR] Opening %s failed: %s\n", data->input_path1, error->message);
      exit (1);
    }
  iter2 = fastq_iter_new (data->input_path2, &error);
  if (error)
    {
      g_printerr ("[ERROR] Opening %s failed: %s\n", data->input_path2, error->message);
      exit (1);
    }

  seq1 = fastq_iter_next (iter1);
  seq2 = fastq_iter_next (iter2);
  while (seq1 != NULL && seq2 != NULL)
    {
      KmerHashNode *node;
      FastqPair    *pair;

      if (!check_sanity (data, seq1, seq2, start, end))
        continue;

      /* Build kmer */
      char_to_bin_prealloc (kmer, seq1->seq + start - 1, k);
      char_to_bin_prealloc (kmer + k_bytes, seq2->seq + start - 1, k);
      node = kmer_hash_table_lookup_or_create (hash_table, kmer);

      /* Insert new pair of fastq sequences */
      if (node->value == NULL)
        {
          pair        = g_slice_new (FastqPair);
          pair->seq1  = fastq_seq_copy (seq1);
          pair->seq2  = fastq_seq_copy (seq2);
          node->value = pair;
        }
      /* Compare and replace */
      else
        {
          pair = (FastqPair*) node->value;
          if (fastq_pair_cmp (data, seq1, seq2, pair->seq1, pair->seq2) > 0)
            {
              fastq_seq_free (pair->seq1);
              fastq_seq_free (pair->seq2);
              pair->seq1 = fastq_seq_copy (seq1);
              pair->seq2 = fastq_seq_copy (seq2);
            }
        }
      seq1 = fastq_iter_next (iter1);
      seq2 = fastq_iter_next (iter2);
    }
  fastq_iter_free (iter1);
  fastq_iter_free (iter2);
  g_free (kmer);

  return hash_table;
}

static KmerHashTable*
remove_clonal (CallbackData   *data,
               KmerHashTable  *hash_table,
               int             round)
{
  KmerHashTableIter iter;
  KmerHashTable    *new_table;
  KmerHashNode     *node;
  unsigned char    *kmer;
  unsigned int      start;
  unsigned int      end;
  unsigned int      k;
  unsigned int      k_bytes;

  if (data->verbose)
    g_printerr ("Filtering clones (round %d)\n", round + 1);

  start     = data->start_coords[round];
  end       = data->end_coords[round];
  k         = end - start + 1;
  k_bytes   = (k + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;
  new_table = kmer_hash_table_new (k_bytes * NUCS_PER_BYTE * 2);
  kmer      = g_malloc0 (2 * k_bytes * sizeof (*kmer));

  kmer_hash_table_iter_init (&iter, hash_table);
  while ((node = kmer_hash_table_iter_next (&iter)) != NULL)
    {
      KmerHashNode *new_node;
      FastqPair    *pair;

      pair = (FastqPair*) node->value;

      if (!check_sanity (data, pair->seq1, pair->seq2, start, end))
        {
          fastq_seq_free (pair->seq1);
          fastq_seq_free (pair->seq2);
          g_slice_free (FastqPair, pair);
          continue;
        }

      /* Build kmer */
      char_to_bin_prealloc (kmer, pair->seq1->seq + start - 1, k);
      char_to_bin_prealloc (kmer + k_bytes, pair->seq2->seq + start - 1, k);
      new_node = kmer_hash_table_lookup_or_create (new_table, kmer);

      /* Insert new pair of fastq sequences */
      if (new_node->value == NULL)
        new_node->value = pair;
      /* Compare and replace */
      else
        {
          FastqPair *new_pair;

          new_pair = (FastqPair*) new_node->value;
          if (fastq_pair_cmp (data, pair->seq1, pair->seq2, new_pair->seq1, new_pair->seq2) > 0)
            {
              fastq_seq_free (new_pair->seq1);
              fastq_seq_free (new_pair->seq2);
              new_pair->seq1 = pair->seq1;
              new_pair->seq2 = pair->seq2;
            }
          else
            {
              fastq_seq_free (pair->seq1);
              fastq_seq_free (pair->seq2);
            }
          g_slice_free (FastqPair, pair);
        }
      node->value = NULL;
    }

  kmer_hash_table_destroy (hash_table);
  g_free (kmer);

  return new_table;
}

static int
fastq_pair_cmp (CallbackData *data,
                FastqSeq     *seq11,
                FastqSeq     *seq12,
                FastqSeq     *seq21,
                FastqSeq     *seq22)
{
  float mean1 = 0;
  float mean2 = 0;
  int   i;

  for (i = 0; i < seq11->size; i++)
    mean1 += seq11->qual[i];
  for (i = 0; i < seq12->size; i++)
    mean1 += seq12->qual[i];
  mean1 /= seq11->size + seq12->size;

  for (i = 0; i < seq21->size; i++)
    mean2 += seq21->qual[i];
  for (i = 0; i < seq22->size; i++)
    mean2 += seq22->qual[i];
  mean2 /= seq21->size + seq22->size;

  if (mean1 == mean2)
    return 0;
  if (mean1 > mean2)
    return 1;
  return -1;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
