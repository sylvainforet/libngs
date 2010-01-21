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
  int                verbose;
  int                print_letter;
  int                print_all;
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
      {"reference", 'r', 0, G_OPTION_ARG_FILENAME,       &data->ref_path,     "Reference genome file", NULL},
      {"bsq",       'b', 0, G_OPTION_ARG_FILENAME_ARRAY, &data->bsq_paths,    "Bsmap bsq file(s)", NULL},
      {"fastq",     'f', 0, G_OPTION_ARG_FILENAME_ARRAY, &data->fastq_paths,  "Fastq file(s)", NULL},
      {"out",       'o', 0, G_OPTION_ARG_FILENAME,       &data->output_path,  "Output file", NULL},
      {"add",       'a', 0, G_OPTION_ARG_FILENAME,       &data->add_path,     "Add results to this file", NULL},
      {"min_qual",  'm', 0, G_OPTION_ARG_INT,            &data->min_qual,     "Minimum base quality", NULL},
      {"verbose",   'v', 0, G_OPTION_ARG_NONE,           &data->verbose,      "Verbose output", NULL},
      {"letter",    'l', 0, G_OPTION_ARG_NONE,           &data->print_letter, "Prepend a column with the letter", NULL},
      {"all",       'w', 0, G_OPTION_ARG_NONE,           &data->print_all,    "Prints all positions (implies l)", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->ref_path     = NULL;
  data->fastq_paths  = NULL;
  data->bsq_paths    = NULL;
  data->output_path  = strdup("-");
  data->add_path     = NULL;
  data->min_qual     = 0;
  data->verbose      = 0;
  data->print_letter = 0;
  data->print_all    = 0;

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

  name_len = strlen (rec->name);
  if (rec->flag == BSQ_MAP_UM &&
      ((rec->name[name_len - 2] == '/' &&
       rec->name[name_len - 1] == '1' &&
       (rec->strand == BSQ_STRAND_W ||
        rec->strand == BSQ_STRAND_C)) ||
      (rec->name[name_len - 2] == '/' &&
       rec->name[name_len - 1] == '2' &&
       (rec->strand == BSQ_STRAND_WC ||
        rec->strand == BSQ_STRAND_CC))))
    {
      SeqDBElement *read_elem;
      SeqDBElement *ref_elem;
      char         *read;
      char         *qual;
      char         *ref;
      unsigned int  start_ref;
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
      read      = data->reads->seqs + read_elem->offset;
      qual      = data->reads->quals + read_elem->offset;
      ref       = data->ref->seqs + start_ref;

      switch (rec->strand)
        {
          case BSQ_STRAND_W:
              /* Dont reverse */
              break;
          case BSQ_STRAND_C:
              /* Reverse ref */
              ref = rev_comp_in_place (ref, read_elem->size);
              is_ref_rev = 1;
              break;
          case BSQ_STRAND_WC:
              /* Reverse read */
              read = rev_comp_in_place (read, read_elem->size);
              qual = rev_in_place (qual, read_elem->size);
              break;
          case BSQ_STRAND_CC:
              /* Reverse both read and ref */
              ref  = rev_comp_in_place (ref, read_elem->size);
              qual = rev_in_place (qual, read_elem->size);
              read = rev_comp_in_place (read, read_elem->size);
              is_ref_rev = 1;
              break;
          default:
              break;
        }

      for (i = 0; i < read_elem->size; i++)
        {
          if (ref[i] == 'C')
            {
              if (qual[i] >= data->min_qual + FASTQ_QUAL_0)
                {
                  unsigned int meth_idx;

                  meth_idx = i;
                  if (is_ref_rev)
                    meth_idx = read_elem->size - i - 1;
                  if (read[i] == 'C')
                    data->counts->meth_index[start_ref + meth_idx]->n_meth++;
                  else if (read[i] == 'T')
                    data->counts->meth_index[start_ref + meth_idx]->n_unmeth++;
                }
            }
        }

      /* Put the reference back into place
       * Should not need to have the reads back in the right orientation ...
       */
      if (is_ref_rev)
        ref  = rev_comp_in_place (ref, read_elem->size);
    }
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
