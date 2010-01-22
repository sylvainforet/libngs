/**
 *
 */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ngs_methylation.h"


typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char              *ref_path;
  char              *input_path;
  char              *output_path;

  SeqDB             *ref;
  RefMethCounts     *counts;

  int                verbose;
  int                min_count;

  unsigned long      n_c;
  unsigned long      n_c_meth;
  unsigned long      n_c_unmeth;
  unsigned long      n_cpg;
  unsigned long      n_cpg_meth;
  unsigned long      n_cpg_unmeth;
  unsigned long      n_chg;
  unsigned long      n_chg_meth;
  unsigned long      n_chg_unmeth;
  unsigned long      n_chh;
  unsigned long      n_chh_meth;
  unsigned long      n_chh_unmeth;
};

static void parse_args    (CallbackData      *data,
                           int               *argc,
                           char            ***argv);

static void load_data     (CallbackData      *data);

static void count_cgs     (CallbackData      *data);

static void cleanup_data  (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;

  parse_args (&data, &argc, &argv);
  load_data (&data);
  count_cgs (&data);
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
      {"reference", 'r', 0, G_OPTION_ARG_FILENAME, &data->ref_path,     "Reference genome file", NULL},
      {"out",       'o', 0, G_OPTION_ARG_FILENAME, &data->output_path,  "Output file", NULL},
      {"verbose",   'v', 0, G_OPTION_ARG_NONE,     &data->verbose,      "Verbose output", NULL},
      {"min_count", 'm', 0, G_OPTION_ARG_INT,      &data->min_count,    "Minimum number of reads", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->ref_path     = NULL;
  data->output_path  = strdup("-");
  data->input_path   = NULL;
  data->min_count    = 0;
  data->verbose      = 0;
  data->n_c          = 0;
  data->n_c_meth     = 0;
  data->n_c_unmeth   = 0;
  data->n_cpg        = 0;
  data->n_cpg_meth   = 0;
  data->n_cpg_unmeth = 0;
  data->n_chg        = 0;
  data->n_chg_meth   = 0;
  data->n_chg_unmeth = 0;
  data->n_chh        = 0;
  data->n_chh_meth   = 0;
  data->n_chh_unmeth = 0;

  context = g_option_context_new ("FILE - Counts various kinds of Cs in the genome");
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
  if (*argc < 2)
    {
      g_printerr ("[ERROR] No input file provided\n");
      exit (1);
    }
  data->input_path = (*argv)[1];
}

static void
load_data (CallbackData *data)
{
  GError *error = NULL;

  data->ref   = seq_db_new ();

  if (data->verbose)
    g_print (">>> Loading reference %s\n", data->ref_path);
  seq_db_load_fasta (data->ref,
                     data->ref_path,
                     &error);
  if (error)
    {
      g_printerr ("[ERROR] Loading reference `%s' failed: %s\n",
                  data->ref_path,
                  error->message);
      exit (1);
    }
  if (data->verbose)
    g_print (">>> Loading meth file %s\n", data->input_path);
  data->counts = ref_meth_counts_load (data->ref,
                                       data->input_path,
                                       &error);
  if (error)
    {
      g_printerr ("[ERROR] Loading meth file `%s' failed: %s\n",
                  data->input_path,
                  error->message);
      exit (1);
    }
}

static void
count_cgs (CallbackData *data)
{
  GHashTableIter  iter;
  GIOChannel     *channel;
  SeqDBElement   *elem;
  char           *buffer;
  GError         *error      = NULL;
  int             use_stdout = 1;

  if (data->verbose)
    g_print (">>> Counting Cs, Gs and CpGs\n");

  g_hash_table_iter_init (&iter, data->ref->index);
  while (g_hash_table_iter_next (&iter, NULL, (gpointer*)&elem))
    {
      guint64       i;
      const guint64 maxi = elem->offset + elem->size - 2;

      /* TODO
       * Take border effects into account?
       * This would mean to also consider:
       * the CpGs that occur at +1 (fromthe start) and -1 (from the end)
       * the Cs that occur at +0 and -0
       */
      for (i = elem->offset + 2; i < maxi; i++)
        {
          if (data->ref->seqs[i] == 'C')
            {
              data->n_c++;
              /* CpG */
              if (data->ref->seqs[i + 1] == 'G')
                {
                  data->n_cpg++;
                  if (data->counts->meth_index[i]->n_meth >= data->min_count)
                    data->n_cpg_meth++;
                  else if (data->counts->meth_index[i]->n_unmeth >= data->min_count)
                    data->n_cpg_unmeth++;
                }
              /* CHG */
              else if (data->ref->seqs[i + 2] == 'G')
                {
                  data->n_chg++;
                  if (data->counts->meth_index[i]->n_meth >= data->min_count)
                    data->n_chg_meth++;
                  else if (data->counts->meth_index[i]->n_unmeth >= data->min_count)
                    data->n_chg_unmeth++;
                }
              /* CHH */
              else
                {
                  data->n_chh++;
                  if (data->counts->meth_index[i]->n_meth >= data->min_count)
                    data->n_chh_meth++;
                  else if (data->counts->meth_index[i]->n_unmeth >= data->min_count)
                    data->n_chh_unmeth++;
                }
            }
          else if (data->ref->seqs[i] == 'G')
            {
              data->n_c++;
              /* CpG */
              if (data->ref->seqs[i - 1] == 'C')
                {
                  data->n_cpg++;
                  if (data->counts->meth_index[i]->n_meth >= data->min_count)
                    data->n_cpg_meth++;
                  else if (data->counts->meth_index[i]->n_unmeth >= data->min_count)
                    data->n_cpg_unmeth++;
                }
              /* CHG */
              else if (data->ref->seqs[i - 2] == 'C')
                {
                  data->n_chg++;
                  if (data->counts->meth_index[i]->n_meth >= data->min_count)
                    data->n_chg_meth++;
                  else if (data->counts->meth_index[i]->n_unmeth >= data->min_count)
                    data->n_chg_unmeth++;
                }
              /* CHH */
              else
                {
                  data->n_chh++;
                  if (data->counts->meth_index[i]->n_meth >= data->min_count)
                    data->n_chh_meth++;
                  else if (data->counts->meth_index[i]->n_unmeth >= data->min_count)
                    data->n_chh_unmeth++;
                }
            }
        }
    }

  if (data->output_path[0] == '-' && data->output_path[1] == '\0')
    channel = g_io_channel_unix_new (STDOUT_FILENO);
  else
    {
      use_stdout = 0;
      channel    = g_io_channel_new_file (data->output_path, "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening output file `%s' failed: %s\n",
                      data->output_path,
                      error->message);
          exit (1);
        }
    }

  buffer = g_strdup_printf ("Methylated and unmethylated Cs and CpGs counts (min count = %d)\n"
                            "  C   total  : %12ld\n"
                            "  CpG total  : %12ld\n"
                            "  CpG meth   : %12ld\n"
                            "  CpG un-meth: %12ld\n"
                            "  CHG total  : %12ld\n"
                            "  CHG meth   : %12ld\n"
                            "  CHG un-meth: %12ld\n"
                            "  CHH total  : %12ld\n"
                            "  CHH meth   : %12ld\n"
                            "  CHH un-meth: %12ld\n",
                            data->min_count,
                            data->n_c,
                            data->n_cpg,
                            data->n_cpg_meth,
                            data->n_cpg_unmeth,
                            data->n_chg,
                            data->n_chg_meth,
                            data->n_chg_unmeth,
                            data->n_chh,
                            data->n_chh_meth,
                            data->n_chh_unmeth);

  g_io_channel_write_chars (channel,
                            buffer,
                            -1,
                            NULL,
                            &error);
  if (error)
    {
      g_printerr ("[ERROR] Writing to output file `%s' failed: %s\n",
                  data->output_path,
                  error->message);
      g_error_free (error);
      error = NULL;
    }
  g_free (buffer);

  if (!use_stdout)
    {
      g_io_channel_shutdown (channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing output file `%s' failed: %s\n",
                      data->output_path,
                      error->message);
          g_error_free (error);
          error = NULL;
        }
    }
  g_io_channel_unref (channel);
}

static void
cleanup_data (CallbackData *data)
{
  if (!data)
    return;
  if (data->ref_path)
    g_free (data->ref_path);
  if (data->output_path)
    g_free (data->output_path);
  if (data->counts)
    ref_meth_counts_destroy (data->counts);
  seq_db_free (data->ref);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
