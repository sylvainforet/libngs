/* Copyright (C) 2010  Sylvain FORET
  data->sanger           = 0;
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
#include <unistd.h>
#include <string.h>

#include "ngs_binseq.h"
#include "ngs_fasta.h"
#include "ngs_fastq.h"
#include "ngs_utils.h"


#define SEED_SIZE 12
#define NB_SEEDS  (1 << (SEED_SIZE * BITS_PER_NUC))


typedef struct _Adaptor Adaptor;

struct _Adaptor
{
  char      *seq;
  char      *rev;
  GPtrArray *names;
  long       counts_fh;
  long       counts_ft;
  long       counts_rh;
  long       counts_rt;
  long       counts;
  int        size;
};

typedef struct _Seed Seed;

struct _Seed
{
  Adaptor     *adaptor;
  Seed        *next;

  unsigned int forward: 1;
  unsigned int head   : 1;
  unsigned int exact  : 1;
};

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char       *adaptors_path;
  char       *reads_path;
  char       *out_path;
  GIOChannel *out_channel;
  GString    *buffer;

  Seed      **seeds;

  long        reads_found;

  Adaptor    *adaptors;
  int         n_adaptors;

  int         len;
  int         mis;
  int         verbose;

  int         use_stdout: 1;
};

static int        iter_func            (FastqSeq       *fastq,
                                        CallbackData   *data);

static void       parse_args           (CallbackData   *data,
                                        int            *argc,
                                        char         ***argv);

static void       load_adaptors        (CallbackData   *data);

static inline int match_head           (CallbackData   *data,
                                        FastqSeq       *fastq,
                                        Adaptor        *adaptor,
                                        int             forward,
                                        int             i,
                                        int             mismatches);

static inline int match_tail           (CallbackData   *data,
                                        FastqSeq       *fastq,
                                        Adaptor        *adaptor,
                                        int             forward,
                                        int             i,
                                        int             mismatches);

static void       init_seeds           (CallbackData   *data);

static void       add_seed             (CallbackData   *data,
                                        Adaptor        *adaptor,
                                        unsigned int    forward,
                                        unsigned int    head);

static guint32    make_seed            (const char     *seq);

static char*      get_adaptor_seed_seq (Adaptor        *adaptor,
                                        unsigned int    forward,
                                        unsigned int    head);

static void       print_summary        (CallbackData   *data);

static void       cleanup_data         (CallbackData   *data);

static int        adaptor_cmp          (const void     *a1,
                                        const void     *a2);

int
main (int    argc,
      char **argv)
{
  CallbackData data;
  GError      *error = NULL;

  parse_args (&data, &argc, &argv);
  load_adaptors (&data);
  init_seeds (&data);
  if (data.verbose)
    g_printerr("Parsing fastq file\n");
  iter_fastq (data.reads_path,
              (FastqIterFunc)iter_func,
              &data,
              &error);
  if (error)
    {
      g_printerr ("[ERROR] Iterating sequences failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
    }
  if (data.verbose)
    print_summary (&data);

  cleanup_data (&data);

  return 0;
}

static void
parse_args (CallbackData   *data,
            int            *argc,
            char         ***argv)
{
  GOptionEntry entries[] =
    {
      {"seqout"  , 'o', 0, G_OPTION_ARG_FILENAME, &data->out_path, "Sequences file"              , NULL},
      {"len",      'l', 0, G_OPTION_ARG_INT,      &data->len,      "Minimum match length"        , NULL},
      {"mis",      'm', 0, G_OPTION_ARG_INT,      &data->mis,      "Maximum number of mismatches", NULL},
      {"verbose",  'v', 0, G_OPTION_ARG_NONE,     &data->verbose,  "Verbose output"              , NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->out_path         = "-";
  data->out_channel      = NULL;
  data->buffer           = NULL;
  data->use_stdout       = 1;
  data->len              = SEED_SIZE;
  data->mis              = 1;
  data->verbose          = 0;
  data->reads_found      = 0;

  context = g_option_context_new ("ADAPTORS READS - Look for adaptors from ADAPTORS (fasta format) in READS (fastq)");
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
      /* TODO look for adaptors in a default location */
      g_printerr ("[ERROR] You must provide an ADAPTORS file and a READS file\n");
      exit (1);
    }
  data->adaptors_path = (*argv)[1];
  data->reads_path    = (*argv)[2];

  if (data->len < SEED_SIZE)
    {
      g_printerr ("[ERROR] Minimum size must be at least: %d\n", SEED_SIZE);
      exit (1);
    }

  /* Set default output to "-" (stdout) */
  if (!data->out_path)
        data->out_path = "-";

  if (data->out_path[0] == '-' && data->out_path[1] == '\0')
    {
      data->use_stdout  = 1;
      data->out_channel = g_io_channel_unix_new (STDOUT_FILENO);
    }
  else
    {
      data->use_stdout  = 0;
      data->out_channel = g_io_channel_new_file (data->out_path, "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening sequence output file failed: %s\n", error->message);
          exit (1);
        }
    }
  data->buffer = g_string_sized_new (512);
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  int     i;
  int     maxi;
  int     ret   = 1;
  int     found = 0;
  char   *mask  = NULL;
  GError *error = NULL;

  if (fastq->size < data->len)
    return ret;

  maxi = fastq->size - SEED_SIZE + 1;

  for (i = 0; i < maxi; i++)
    {
      Seed         *seed;
      const guint32 idx = make_seed (fastq->seq + i);

      if (data->seeds[idx] == NULL)
        continue;

      for (seed = data->seeds[idx]; seed != NULL; seed = seed->next)
        {
          int mismatches = seed->exact ? 0 : 1;
          int match      = 0;

          if (seed->head)
            {
              match = match_head (data,
                                  fastq,
                                  seed->adaptor,
                                  seed->forward,
                                  i,
                                  mismatches);
              if (match)
                {
                  int j;
                  int maxj;

                  if (seed->forward)
                    seed->adaptor->counts_fh++;
                  else
                    seed->adaptor->counts_rh++;

                  if (mask == NULL)
                    mask = g_malloc0 (fastq->size);

                  maxj = fastq->size;
                  if (i + seed->adaptor->size < fastq->size)
                    maxj = i + seed->adaptor->size;
                  for (j = i; j < maxj; j++)
                    mask[j] = 1;
                }
            }
          else
            {
              match = match_tail (data,
                                  fastq,
                                  seed->adaptor,
                                  seed->forward,
                                  i,
                                  mismatches);
              if (match)
                {
                  int j;
                  int minj;

                  if (seed->forward)
                    seed->adaptor->counts_ft++;
                  else
                    seed->adaptor->counts_rt++;

                  if (mask == NULL)
                    mask = g_malloc0 (fastq->size);

                  minj = 0;
                  if (i - seed->adaptor->size + 1 > 0)
                    minj = i - seed->adaptor->size + 1;
                  for (j = minj; j <= i; j++)
                    mask[j] = 1;
                }
            }
          if (match)
            found = 1;
        }
    }
  if (found)
    {
      int current_start = -1;
      int current_len   = 0;;
      int max_start     = -1;
      int max_len       = 0;

      for (i = 0; i < fastq->size; i++)
        {
          if (mask[i] == 0)
            {
              if (current_start == -1)
                  current_start = i;
              current_len++;
            }
          else
            {
              if (current_start != -1)
                {
                  if (current_len > max_len)
                    {
                      max_start = current_start;
                      max_len   = current_len;
                    }
                  current_start = -1;
                  current_len   = 0;
                }
            }
        }
      if (current_start != -1 && current_len > max_len)
        {
          max_start = current_start;
          max_len   = current_len;
        }

      if (max_start == -1)
        fastq_write (data->out_channel,
                     data->buffer,
                     fastq->name,
                     "",
                     "",
                   &error);
      else
        fastq_write_fragment (data->out_channel,
                              data->buffer,
                              fastq->name,
                              fastq->seq,
                              fastq->qual,
                              max_start,
                              max_start + max_len,
                              &error);
      data->reads_found++;
      g_free (mask);

      if (error != NULL)
        {
          g_printerr ("[ERROR] Writing sequence failed: %s\n",
                      error->message);
          g_error_free (error);
          return 0;
        }
    }
  else
    {
      fastq_write (data->out_channel,
                   NULL,
                   fastq->name,
                   fastq->seq,
                   fastq->qual,
                   &error);
      if (error != NULL)
        {
          g_printerr ("[ERROR] Writing sequence failed: %s\n",
                      error->message);
          g_error_free (error);
          return 0;
        }
    }

  return ret;
}

static inline int
match_head (CallbackData *data,
            FastqSeq     *fastq,
            Adaptor      *adaptor,
            int           forward,
            int           i,
            int           mismatches)
{
  if (i + data->len > fastq->size)
    return 0;
  else
    {
      const char *seq = forward ? adaptor->seq : adaptor->rev;
      int         j   = SEED_SIZE;
      int         k   = SEED_SIZE + i;

      while (j < adaptor->size && k < fastq->size)
        {
          if (seq[j] != fastq->seq[k])
            {
              ++mismatches;
              if (mismatches > data->mis)
                return 0;
            }
          ++j;
          ++k;
        }
      if (mismatches <= data->mis)
        {
          /* g_printerr ("H %d: %s %s (%d - %d) / (%d - %d)\n", forward, fastq->seq, adaptor->seq, i, k - 1, 0, j - 1); */
          return 1;
        }
    }

  return 0;
}

static inline int
match_tail (CallbackData *data,
            FastqSeq     *fastq,
            Adaptor      *adaptor,
            int           forward,
            int           i,
            int           mismatches)
{
  if (i + SEED_SIZE - data->len < 0)
    return 0;
  else
    {
      const char *seq = forward ? adaptor->seq : adaptor->rev;
      int         j   = adaptor->size - SEED_SIZE - 1;
      int         k   = i - 1;

      while (j >= 0 && k >= 0)
        {
          if (seq[j] != fastq->seq[k])
            {
              ++mismatches;
              if (mismatches > data->mis)
                return 0;
            }
          --j;
          --k;
        }
      /* Full match: counted as head */
      if (j == -1)
        return 0;
      if (mismatches <= data->mis)
        {
          /* g_printerr ("T %d: %s %s (%d - %d) / (%d - %d)\n", forward, fastq->seq, adaptor->seq, k + 1, i + SEED_SIZE - 1, j + 1, adaptor->size - 1); */
          return 1;
        }
    }

  return 0;
}

static void
load_adaptors (CallbackData *data)
{
  GHashTableIter iter_ht;
  GHashTable    *adapt_ht;
  FastaIter     *iter;
  FastaSeq      *seq;
  char          *adapt_key;
  GPtrArray     *adapt_val;
  GError        *error = NULL;
  int            i = 0;

  if (data->verbose)
    g_printerr ("Loading adaptors\n");

  iter = fasta_iter_new (data->adaptors_path, &error);

  if (error)
    {
      g_printerr ("[ERROR] Iterating adaptor sequences failed: %s\n", error->message);
      exit (1);
    }

  adapt_ht = g_hash_table_new (g_str_hash, g_str_equal);

  for (seq = fasta_iter_next (iter) ; seq != NULL; seq = fasta_iter_next (iter))
    {
      GPtrArray *names;

      if (g_hash_table_lookup_extended (adapt_ht, seq->seq, NULL, (gpointer*)&names))
        g_ptr_array_add (names, g_strdup (seq->name));
      else
        {
          rev_comp_in_place (seq->seq, seq->size);

          if (g_hash_table_lookup_extended (adapt_ht, seq->seq, NULL, (gpointer*)&names))
            g_ptr_array_add (names, g_strdup_printf ("%s (-)", seq->name));
          else
            {
              rev_comp_in_place (seq->seq, seq->size);
              names = g_ptr_array_new_with_free_func (g_free);
              g_ptr_array_add (names, g_strdup (seq->name));
              g_hash_table_insert (adapt_ht, g_strdup (seq->seq), names);
            }
        }
    }

  data->n_adaptors = g_hash_table_size (adapt_ht);
  data->adaptors   = g_malloc (data->n_adaptors * sizeof (*data->adaptors));

  g_hash_table_iter_init (&iter_ht, adapt_ht);
  while (g_hash_table_iter_next (&iter_ht, (gpointer*)&adapt_key, (gpointer*)&adapt_val)) 
    {
      data->adaptors[i].seq       = adapt_key;
      data->adaptors[i].rev       = g_strdup(adapt_key);
      data->adaptors[i].names     = adapt_val;
      data->adaptors[i].counts_fh = 0;
      data->adaptors[i].counts_ft = 0;
      data->adaptors[i].counts_rh = 0;
      data->adaptors[i].counts_rt = 0;
      data->adaptors[i].counts    = 0;
      data->adaptors[i].size      = strlen(adapt_key);
      rev_comp_in_place (data->adaptors[i].rev,
                         data->adaptors[i].size);
      ++i;
    }
  g_hash_table_destroy (adapt_ht);
}

static void
init_seeds (CallbackData *data)
{
  int i;

  if (data->verbose)
    g_printerr ("Initialiasing seeds\n");

  data->seeds = g_malloc0 (NB_SEEDS * sizeof (*data->seeds));

  for (i = 0; i < data->n_adaptors; i++)
    {
      if (data->adaptors[i].size < data->len)
        {
          if (data->verbose)
            g_printerr ("[WARNING] Adaptor %s too small, ignoring\n", data->adaptors[i].seq);

          continue;
        }

      add_seed (data, data->adaptors + i, 0, 0);
      add_seed (data, data->adaptors + i, 1, 0);
      add_seed (data, data->adaptors + i, 0, 1);
      add_seed (data, data->adaptors + i, 1, 1);
    }
}

static void
add_seed (CallbackData *data,
          Adaptor      *adaptor,
          unsigned int  forward,
          unsigned int  head)
{
  Seed    *seed;
  char    *seq;
  guint32  idx;
  int      i;

  seq = get_adaptor_seed_seq (adaptor, forward, head);
  idx = make_seed (seq);

  seed                    = g_slice_new0 (Seed);
  seed->adaptor           = adaptor;
  seed->forward           = forward;
  seed->head              = head;
  seed->exact             = 1;
  seed->next              = data->seeds[idx];
  data->seeds[idx]        = seed;

  for (i = 0; i < SEED_SIZE; i++)
    {
      const char c = seq[i];
      int j;

      for (j = 0; j < 4; j++)
        {
          const char d = bin_to_char_table[j];

          if (d == c)
            continue;

          seq[i]                        = d;
          idx                           = make_seed (seq);
          seed                          = g_slice_new0 (Seed);
          seed->adaptor                 = adaptor;
          seed->forward                 = forward;
          seed->head                    = head;
          seed->exact                   = 0;
          seed->next                    = data->seeds[idx];
          data->seeds[idx]              = seed;
        }
      seq[i] = c;
    }
}

static guint32
make_seed (const char *seq)
{
  guint32 seed = 0;
  int     i;

  for (i = 0; i < SEED_SIZE; i++)
    {
      seed <<= 2;
      seed |= char_to_bin_table[(int)seq[i]];
    }
  return seed;
}

static char*
get_adaptor_seed_seq (Adaptor     *adaptor,
                      unsigned int forward,
                      unsigned int head)
{
  char *seq;

  seq = adaptor->seq;
  if (!forward)
    seq = adaptor->rev;

  if (head)
    return seq;

  return seq + adaptor->size - SEED_SIZE;
}


static void
print_summary (CallbackData *data)
{
  Adaptor **adaptors;
  int       i;

  adaptors = g_malloc (data->n_adaptors * sizeof (*adaptors));
  for (i = 0; i < data->n_adaptors; i++)
    {
      adaptors[i]         = data->adaptors + i;
      adaptors[i]->counts = adaptors[i]->counts_fh +
                            adaptors[i]->counts_ft +
                            adaptors[i]->counts_rh +
                            adaptors[i]->counts_rt;
    }

  qsort (adaptors,
         data->n_adaptors,
         sizeof (*adaptors),
         adaptor_cmp);

  g_printerr ("*** Summary of adaptors found\n");
  for (i = 0; i < data->n_adaptors; i++)
    {
      if (adaptors[i]->counts > 0)
        {
          int j;

          g_printerr ("* Adaptor %s", adaptors[i]->seq);

          for (j = 0; j < adaptors[i]->names->len; j++)
            g_printerr (" %s", (char*)g_ptr_array_index (adaptors[i]->names, j));

          g_printerr ("\n");
          g_printerr ("Found %ld times (fh:%ld, ft:%ld, rh:%ld, rt:%ld)\n",
                      adaptors[i]->counts,
                      adaptors[i]->counts_fh,
                      adaptors[i]->counts_ft,
                      adaptors[i]->counts_rh,
                      adaptors[i]->counts_rt);
        }
    }
  g_printerr ("*** Reads found with adaptors: %ld\n", data->reads_found);
}

static void
cleanup_data (CallbackData *data)
{
  int i;

  if (data->verbose)
    g_printerr ("Cleaning up\n");

  if (data->out_channel)
    {
      if (!data->use_stdout)
        {
          GError *error = NULL;

          g_io_channel_shutdown (data->out_channel, TRUE, &error);
          if (error)
            {
              g_printerr ("[ERROR] Closing sequence output file failed: %s\n", error->message);
              g_error_free (error);
              error = NULL;
            }
        }
      g_io_channel_unref (data->out_channel);
    }

  if (data->adaptors != NULL)
    {
      for (i = 0; i < data->n_adaptors; i++)
        {
          g_ptr_array_free (data->adaptors[i].names, TRUE);
          g_free (data->adaptors[i].seq);
          g_free (data->adaptors[i].rev);
        }

      g_free (data->adaptors);
      data->adaptors   = NULL;
      data->n_adaptors = 0;
    }

  if (data->seeds != NULL)
    {
      for (i = NB_SEEDS - 1; i >= 0; --i)
        {
          if (data->seeds[i] != NULL)
            {
              Seed *tmp;

              for (tmp = data->seeds[i]; tmp != NULL; )
                {
                  tmp = data->seeds[i]->next;
                  g_slice_free (Seed, data->seeds[i]);
                  data->seeds[i] = tmp;
                }
            }
        }
      g_free (data->seeds);
      data->seeds = NULL;
    }

  if (data->buffer != NULL)
    {
      g_string_free (data->buffer, TRUE);
      data->buffer = NULL;
    }
}

static int
adaptor_cmp (const void     *a1,
             const void     *a2)
{
  Adaptor *adaptor1 = *(Adaptor**)a1;
  Adaptor *adaptor2 = *(Adaptor**)a2;

  if (adaptor1->counts < adaptor2->counts)
    return 1;
  if (adaptor1->counts > adaptor2->counts)
    return -1;
  return 0;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
