/**
 *
 */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ngs_methylation.h"

typedef enum
{
  METH_CPG = 0,
  METH_CHG,
  METH_CHH,
  METH_ALL,
  METH_TYPE_NB
}
MethylationType;

char *meth_type_names[METH_TYPE_NB] =
{
  [METH_CPG] = "cpg",
  [METH_CHG] = "chg",
  [METH_CHH] = "chh",
  [METH_ALL] = "all"
};


typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char              *ref_path;
  char              *input_path;
  char              *output_path;
  char              *meth_type_str;

  SeqDB             *ref;
  RefMethCounts     *counts;

  int                verbose;
  int                ratio;
  int                min_count_meth;
  int                min_count_unmeth;
  int                min_count_tot;
  int                meth_type;
};

static void parse_args    (CallbackData      *data,
                           int               *argc,
                           char            ***argv);

static void load_data     (CallbackData      *data);

static void write_ratios  (CallbackData      *data);

static void cleanup_data  (CallbackData      *data);

static int  get_meth_type (const char        *str);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;

  parse_args (&data, &argc, &argv);
  load_data (&data);
  write_ratios (&data);
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
      {"reference",        'r', 0, G_OPTION_ARG_FILENAME, &data->ref_path,         "Reference genome file", NULL},
      {"out",              'o', 0, G_OPTION_ARG_FILENAME, &data->output_path,      "Output file", NULL},
      {"verbose",          'v', 0, G_OPTION_ARG_NONE,     &data->verbose,          "Verbose output", NULL},
      {"ratio",            'a', 0, G_OPTION_ARG_NONE,     &data->ratio,            "Output ratio instead of meth/unmeth numbers (default)", NULL},
      {"min_count_meth",   'm', 0, G_OPTION_ARG_INT,      &data->min_count_meth,   "Minimum number of methylated reads", NULL},
      {"min_count_unmeth", 'u', 0, G_OPTION_ARG_INT,      &data->min_count_unmeth, "Minimum number of unmethylated reads", NULL},
      {"min_count_tot",    't', 0, G_OPTION_ARG_INT,      &data->min_count_tot,    "Minimum total number of reads", NULL},
      {"meth_type",        'c', 0, G_OPTION_ARG_STRING,   &data->meth_type_str,    "Type of methylation", "[cpg|chg|chh|all]"},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->ref_path         = NULL;
  data->output_path      = strdup("-");
  data->input_path       = NULL;
  data->min_count_meth   = 0;
  data->min_count_unmeth = 0;
  data->min_count_tot    = 0;
  data->verbose          = 0;
  data->ratio            = 0;
  data->meth_type_str    = NULL;

  context = g_option_context_new ("FILE - Prints the methylation ratios");
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
  if (!data->meth_type_str)
    data->meth_type_str = strdup (meth_type_names[METH_CPG]);
  data->meth_type = get_meth_type (data->meth_type_str);
  if (data->meth_type < 0)
    {
      g_printerr ("[ERROR] meth_type (-c) must be one of: %s, %s, %s, %s\n",
                  meth_type_names[METH_CPG],
                  meth_type_names[METH_CHG],
                  meth_type_names[METH_CHH],
                  meth_type_names[METH_ALL]);
      exit (1);
    }
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
write_ratios (CallbackData *data)
{
  GHashTableIter  iter;
  GIOChannel     *channel;
  SeqDBElement   *elem;
  GString        *buffer;
  GError         *error      = NULL;
  int             use_stdout = 1;

  if (data->verbose)
    g_print (">>> Writing ratios\n");

  buffer = g_string_new (NULL);
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
          int print_this_one = 0;

          if (data->ref->seqs[i] == 'C')
            {
              if (data->meth_type == METH_ALL)
                print_this_one = 1;
              /* CpG */
              else if (data->ref->seqs[i + 1] == 'G')
                {
                  if (data->meth_type == METH_CPG)
                    print_this_one = 1;
                }
              /* CHG */
              else if (data->ref->seqs[i + 2] == 'G')
                {
                  if (data->meth_type == METH_CHG)
                    print_this_one = 1;
                }
              /* CHH */
              else
                {
                  if (data->meth_type == METH_CHH)
                    print_this_one = 1;
                }
            }
          else if (data->ref->seqs[i] == 'G')
            {
              if (data->meth_type == METH_ALL)
                print_this_one = 1;
              /* CpG */
              else if (data->ref->seqs[i - 1] == 'C')
                {
                  if (data->meth_type == METH_CPG)
                    print_this_one = 1;
                }
              /* CHG */
              else if (data->ref->seqs[i - 2] == 'C')
                {
                  if (data->meth_type == METH_CHG)
                    print_this_one = 1;
                }
              /* CHH */
              else
                {
                  if (data->meth_type == METH_CHH)
                    print_this_one = 1;
                }
            }
          if (print_this_one &&
              data->counts->meth_index[i]->n_meth >= data->min_count_meth &&
              data->counts->meth_index[i]->n_unmeth >= data->min_count_unmeth &&
              data->counts->meth_index[i]->n_meth + data->counts->meth_index[i]->n_unmeth >= data->min_count_tot)
            {
              if (data->ratio)
                g_string_append_printf (buffer, "%.3f\n",
                                        ((float)data->counts->meth_index[i]->n_meth) /
                                        (data->counts->meth_index[i]->n_meth + data->counts->meth_index[i]->n_unmeth));
              else
                g_string_append_printf (buffer, "%d\t%d\n",
                                        data->counts->meth_index[i]->n_meth,
                                        data->counts->meth_index[i]->n_unmeth);
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
  g_io_channel_write_chars (channel,
                            buffer->str,
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
  g_string_free (buffer, TRUE);

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
  if (data->meth_type_str)
    g_free (data->meth_type_str);
  seq_db_free (data->ref);
}

static int
get_meth_type (const char *str)
{
  int i;

  if (!str || !*str)
    return -1;

  for (i = 0; i < METH_TYPE_NB; i++)
    if (!strcmp (str, meth_type_names[i]))
      return i;

  return -1;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
