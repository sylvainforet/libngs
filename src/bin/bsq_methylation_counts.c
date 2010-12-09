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

#include "ngs_bsq.h"
#include "ngs_fasta.h"
#include "ngs_fastq.h"
#include "ngs_methylation.h"
#include "ngs_seq_db.h"
#include "ngs_utils.h"


typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char              *ref_path;
  char             **fastq_paths;
  char             **bsq_paths;
  char              *output_path;
  char              *add_path;

  SeqDB             *reads;
  SeqDB             *ref;
  RefMethCounts     *counts;

  int                min_qual;
  int                max_chh;
  int                verbose;
  int                check_strand;
  int                print_letter;
  int                print_all;
  int                sanger;

  unsigned long int  n_chh_filtered;
  unsigned long int  n_bad_orientation;
  unsigned long int  n_nqs_filtered;

  unsigned int       trim_tag;
  unsigned int       nqs_each_side;
  unsigned int       nqs_mismatches;
  unsigned int       nqs_neighbor_qual;
};

static void parse_args    (CallbackData      *data,
                           int               *argc,
                           char            ***argv);

static void load_data     (CallbackData      *data);

static void map_data      (CallbackData      *data);

static int  iter_bsq_func (BsqRecord         *rec,
                           CallbackData      *data);

static void cleanup_data  (CallbackData      *data);


int
main (int    argc,
      char **argv)
{
  CallbackData  data;
  GError       *error = NULL;

  parse_args (&data, &argc, &argv);
  if (data.verbose)
    g_print (">>> Loading Data\n");
  load_data (&data);
  if (data.add_path)
    {
      ref_meth_counts_add_path (data.counts,
                                data.ref,
                                data.add_path,
                                &error);
      if (error)
        {
          g_printerr ("[ERROR] failed to load meth file `%s': %s\n",
                      data.add_path,
                      error->message);
          exit (1);
        }
    }
  if (data.verbose)
    g_print (">>> Mapping Data\n");
  map_data (&data);
  if (data.verbose)
    g_print (">>> Writing Data\n");
  ref_meth_counts_write (data.counts,
                         data.ref,
                         data.output_path,
                         data.print_letter,
                         data.print_all,
                         &error);
  if (error)
    {
      g_printerr ("[ERROR] failed to write meth file `%s': %s\n",
                  data.output_path,
                  error->message);
      exit (1);
    }
  if (data.verbose)
    {
      g_print ("Number of CHH filtered reads    : %lu\n", data.n_chh_filtered);
      g_print ("Number of bad orientation reads : %lu\n", data.n_bad_orientation);
      g_print ("Number of NQS filtered bases    : %lu\n", data.n_nqs_filtered);
    }
  cleanup_data (&data);

  return 0;
}

static void
parse_args (CallbackData      *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      /* Input options */
      {"reference", 'r', 0, G_OPTION_ARG_FILENAME,       &data->ref_path,    "Reference genome file", NULL},
      {"bsq",       'b', 0, G_OPTION_ARG_FILENAME_ARRAY, &data->bsq_paths,   "Bsmap bsq file(s)", NULL},
      {"fastq",     'f', 0, G_OPTION_ARG_FILENAME_ARRAY, &data->fastq_paths, "Fastq file(s)", NULL},
      {"out",       'o', 0, G_OPTION_ARG_FILENAME,       &data->output_path, "Output file", NULL},
      {"add",       'a', 0, G_OPTION_ARG_FILENAME,       &data->add_path,    "Add results to this file", NULL},
      {"sanger",    'g', 0, G_OPTION_ARG_NONE,           &data->sanger,      "Qualities in sanger format", NULL},

      /* Mapping options */
      {"min_qual",     'm', 0, G_OPTION_ARG_INT,  &data->min_qual,     "Minimum base quality", NULL},
      {"max_chh",      'M', 0, G_OPTION_ARG_INT,  &data->max_chh,      "Maximum number of consecutive CHH per read", NULL},
      {"trim_tag",     'T', 0, G_OPTION_ARG_INT,  &data->trim_tag,     "Number of tag letters to trim", NULL},
      {"check_strand", 'c', 0, G_OPTION_ARG_NONE, &data->check_strand, "Check strand orientation", NULL},

      {"nqs_each_side",     's', 0, G_OPTION_ARG_INT, &data->nqs_each_side,     "Number of bases to consider on each side for the NQS criterion", NULL},
      {"nqs_mismatches",    't', 0, G_OPTION_ARG_INT, &data->nqs_mismatches,    "Maximum number of mismatches for the NQS criterion", NULL},
      {"nqs_neighbor_qual", 'n', 0, G_OPTION_ARG_INT, &data->nqs_neighbor_qual, "Minimum neighbor quality for the NQS criterion", NULL},

      /* Output options */
      {"verbose", 'v', 0, G_OPTION_ARG_NONE, &data->verbose,      "Verbose output", NULL},
      {"letter",  'l', 0, G_OPTION_ARG_NONE, &data->print_letter, "Prepend a column with the letter", NULL},
      {"all",     'w', 0, G_OPTION_ARG_NONE, &data->print_all,    "Prints all positions (implies l)", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->ref_path          = NULL;
  data->fastq_paths       = NULL;
  data->bsq_paths         = NULL;
  data->output_path       = strdup("-");
  data->add_path          = NULL;
  data->sanger            = 0;
  data->min_qual          = 0;
  data->max_chh           = 0;
  data->check_strand      = 0;
  data->verbose           = 0;
  data->print_letter      = 0;
  data->print_all         = 0;
  data->trim_tag          = 0;
  data->n_chh_filtered    = 0;
  data->n_bad_orientation = 0;
  data->n_nqs_filtered    = 0;
  data->nqs_each_side     = 0;
  data->nqs_mismatches    = 0;
  data->nqs_neighbor_qual = 0;

  context = g_option_context_new (" - Finds all the potentially methylated Cs");
  g_option_context_add_group (context, get_fasta_option_group ());
  g_option_context_add_group (context, get_fastq_option_group ());
  g_option_context_add_group (context, get_bsq_option_group ());
  g_option_context_add_main_entries (context, entries, NULL);
  if (!g_option_context_parse (context, argc, argv, &error))
    {
      g_printerr ("[ERROR] Option parsing failed: %s\n", error->message);
      exit (1);
    }
  g_option_context_free (context);

  if (!data->ref_path)
    {
      g_printerr ("[ERROR] You must specify a reference genome with -r\n");
      exit (1);
    }
  if (!data->fastq_paths || !*data->fastq_paths)
    {
      g_printerr ("[ERROR] You must specify one or more fastq files with -f\n");
      exit (1);
    }
  if (!data->bsq_paths || !*data->fastq_paths)
    {
      g_printerr ("[ERROR] You must specify one or more bsq files with -b\n");
      exit (1);
    }
  if (data->min_qual >= 35)
    g_printerr ("[WARNING] Minimum quality is set to %d, "
                "it is likely that few or no bases will be considered\n",
                data->min_qual);

  if (data->sanger)
    {
      data->min_qual          += FASTQ_QUAL_0_SANGER;
      data->nqs_neighbor_qual += FASTQ_QUAL_0_SANGER;
    }
  else
    {
      data->min_qual          += FASTQ_QUAL_0;
      data->nqs_neighbor_qual += FASTQ_QUAL_0;
    }

  if (data->print_all)
    data->print_letter = 1;
}

static void
load_data (CallbackData *data)
{
  GError         *error = NULL;
  char          **tmp;

  data->ref   = seq_db_new ();
  data->reads = seq_db_new ();

  if (data->verbose)
    g_print (">>> Loading Fasta\n");
  seq_db_load_fasta (data->ref, data->ref_path, &error);
  if (error)
    {
      g_printerr ("[ERROR] Loading reference failed: %s\n", error->message);
      exit (1);
    }
  if (data->verbose)
    g_print (">>> Loading Fastq\n");
  for (tmp = data->fastq_paths; *tmp; tmp++)
    {
      seq_db_load_fastq (data->reads, *tmp, &error);
      if (error)
        {
          g_printerr ("[ERROR] Loading reads failed: %s\n", error->message);
          exit (1);
        }
    }
  data->counts = ref_meth_counts_create (data->ref);
}

static int
iter_bsq_func (BsqRecord    *rec,
               CallbackData *data)
{
  int name_len;
  int strand_ok = 1;

  name_len = strlen (rec->name);
  if (data->check_strand)
    {
      strand_ok = 0;
      if (rec->flag == BSQ_MAP_UM &&
          ((rec->name[name_len - 2] == '/' &&
           rec->name[name_len - 1] == '1' &&
           (rec->strand == BSQ_STRAND_W ||
            rec->strand == BSQ_STRAND_C)) ||
          (rec->name[name_len - 2] == '/' &&
           rec->name[name_len - 1] == '2' &&
           (rec->strand == BSQ_STRAND_WC ||
            rec->strand == BSQ_STRAND_CC))))
        strand_ok = 1;
    }
  if (strand_ok)
    {
      SeqDBElement *read_elem;
      SeqDBElement *ref_elem;
      char         *read;
      char         *qual;
      char         *ref;
      unsigned int  start_ref;
      unsigned int  read_size;
      unsigned int  i;
      int           is_ref_rev = 0;

      read_elem = g_hash_table_lookup (data->reads->index, rec->name);
      if (!read_elem)
        {
          g_printerr ("[WARNING] Read `%s' not found\n", rec->name);
          return 1;
        }
      ref_elem = g_hash_table_lookup (data->ref->index, rec->ref);
      if (!ref_elem)
        {
          g_printerr ("[WARNING] Reference `%s' not found\n", rec->ref);
          return 1;
        }

      start_ref = ref_elem->offset + rec->loc - 1;
      read      = data->reads->seqs + read_elem->offset + data->trim_tag;
      qual      = data->reads->quals + read_elem->offset + data->trim_tag;
      ref       = data->ref->seqs + start_ref;
      read_size = read_elem->size - data->trim_tag;

      switch (rec->strand)
        {
          case BSQ_STRAND_W:
              /* Dont reverse */
              break;
          case BSQ_STRAND_C:
              /* Reverse ref */
              ref = rev_comp_in_place (ref, read_size);
              is_ref_rev = 1;
              break;
          case BSQ_STRAND_WC:
              /* Reverse read */
              read = rev_comp_in_place (read, read_size);
              qual = rev_in_place (qual, read_size);
              break;
          case BSQ_STRAND_CC:
              /* Reverse both read and ref */
              ref  = rev_comp_in_place (ref, read_size);
              qual = rev_in_place (qual, read_size);
              read = rev_comp_in_place (read, read_size);
              is_ref_rev = 1;
              break;
          default:
              break;
        }

      /* Dont take into account reads that have more than three consecutive CHH
       * See Cokus et al, Nature, 2008.
       */
      if (data->max_chh > 0)
        {
          const unsigned int maxi = read_size - 2;
          unsigned int       j;
          unsigned int       n_chh = 0;

          for (i = 0; i < maxi; i++)
            {
              const unsigned int maxj = MIN (i + 9, read_size);

              n_chh = 0;
              for (j = i; j < maxj; j++)
                {
                  if (read[j] == 'C' && qual[j] >= data->min_qual)
                    {
                      if (read[j + 1] != 'G' && read[j + 2] != 'G')
                        n_chh++;
                      else
                        n_chh = 0;
                    }
                  if (n_chh == data->max_chh)
                    {
                      if (data->verbose)
                        g_print ("CHH-filtered %s on %s\n",
                                 read_elem->name,
                                 ref_elem->name);
                      data->n_chh_filtered++;
                      goto reverse;
                    }
                }
            }
        }

      /* No NQS */
      if (data->nqs_each_side <= 0)
        for (i = 0; i < read_size; i++)
          {
            if (ref[i] == 'C')
              {
                if (qual[i] >= data->min_qual)
                  {
                    unsigned int meth_idx;

                    meth_idx = i;
                    if (is_ref_rev)
                      meth_idx = read_size - i - 1;
                    if (read[i] == 'C')
                      data->counts->meth_index[start_ref + meth_idx]->n_meth++;
                    else if (read[i] == 'T')
                      data->counts->meth_index[start_ref + meth_idx]->n_unmeth++;
                  }
              }
          }
      /* NQS */
      else
        {
          const unsigned int start = data->nqs_each_side;
          const unsigned int end   = read_size - data->nqs_each_side;

          for (i = start; i < end; i++)
            {
              if (ref[i] == 'C')
                {
                  unsigned int j;
                  unsigned int mismatches = 0;
                  int          good       = 1;

                  for (j = 1; j <= data->nqs_each_side; j++)
                    {
                      /* Quality */
                      if (qual[i - j] < data->nqs_neighbor_qual ||
                          qual[i + j] < data->nqs_neighbor_qual)
                        {
                          good = 0;
                          break;
                        }
                      if (ref[i - j] == 'C')
                        {
                          if (read[i - j] != 'C' && read[i - j] != 'T')
                            {
                              ++mismatches;
                              if (mismatches > data->nqs_mismatches)
                                {
                                  good = 0;
                                  break;
                                }
                            }
                        }
                      else if (ref[i - j] != read[i - j])
                        {
                          ++mismatches;
                          if (mismatches > data->nqs_mismatches)
                            {
                              good = 0;
                              break;
                            }
                        }
                      if (ref[i + j] == 'C')
                        {
                          if (read[i + j] != 'C' && read[i + j] != 'T')
                            {
                              ++mismatches;
                              if (mismatches > data->nqs_mismatches)
                                {
                                  good = 0;
                                  break;
                                }
                            }
                        }
                      else if (ref[i + j] != read[i + j])
                        {
                          ++mismatches;
                          if (mismatches > data->nqs_mismatches)
                            {
                              good = 0;
                              break;
                            }
                        }
                    }
                  if (!good)
                    {
                      data->n_nqs_filtered++;
                      continue;
                    }

                  if (qual[i] >= data->min_qual)
                    {
                      unsigned int meth_idx;

                      meth_idx = i;
                      if (is_ref_rev)
                        meth_idx = read_size - i - 1;
                      if (read[i] == 'C')
                        data->counts->meth_index[start_ref + meth_idx]->n_meth++;
                      else if (read[i] == 'T')
                        data->counts->meth_index[start_ref + meth_idx]->n_unmeth++;
                    }
                }
            }
        }

      /* Put the reference back into place
       * Should not need to have the reads back in the right orientation ...
       */
reverse:
      if (is_ref_rev)
        ref  = rev_comp_in_place (ref, read_size);
    }
  else
    data->n_bad_orientation++;

  return 1;
}

static void
map_data (CallbackData *data)
{
  char  **tmp;

  for (tmp = data->bsq_paths; *tmp; tmp++)
    {
      GError *error = NULL;

      iter_bsq (*tmp,
                (BsqIterFunc)iter_bsq_func,
                data,
                &error);
      if (error)
        {
          g_printerr ("[ERROR] Loading bsq file `%s' failed: %s\n",
                      *tmp, error->message);
          exit (1);
        }
    }
}

static void
cleanup_data (CallbackData *data)
{
  if (!data)
    return;
  if (data->fastq_paths)
    g_strfreev (data->fastq_paths);
  if (data->bsq_paths)
    g_strfreev (data->bsq_paths);
  if (data->ref_path)
    g_free (data->ref_path);
  if (data->output_path)
    g_free (data->output_path);
  if (data->counts)
    g_free (data->counts);
  seq_db_free (data->ref);
  seq_db_free (data->reads);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
