/**
 *
 */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ngs_seq_db.h"

typedef struct _MethData MethData;

struct _MethData
{
  unsigned int n_meth;
  unsigned int n_unmeth;

  /*
  unsigned int n_quals
  char        *quals;
  */
};

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char              *ref_path;
  char              *output_path;

  SeqDB             *ref;
  MethData         **meth_index;
  MethData          *meth_data;

  int                verbose;
  int                print_letter;
  int                print_all;
};

static void parse_args    (CallbackData      *data,
                           int               *argc,
                           char            ***argv);

static void load_ref      (CallbackData      *data);

static void load_previous (CallbackData      *data,
                           const char        *path);

static void write_meth    (CallbackData      *data);

static void cleanup_data  (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;
  int           i;

  parse_args (&data, &argc, &argv);

  load_ref (&data);
  for (i = 1; i < argc; i++)
    load_previous (&data, argv[i]);
  write_meth (&data);
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
      {"letter",    'l', 0, G_OPTION_ARG_NONE,     &data->print_letter, "Prepend a column with the letter", NULL},
      {"all",       'w', 0, G_OPTION_ARG_NONE,     &data->print_all,    "Prints all positions (implies l)", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->ref_path     = NULL;
  data->output_path  = strdup("-");
  data->verbose      = 0;
  data->print_letter = 0;
  data->print_all    = 0;

  context = g_option_context_new ("FILE ... - Merges a set of CG files");
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

  if (data->print_all)
    data->print_letter = 1;
}

static void
load_ref (CallbackData *data)
{
  GError         *error = NULL;
  unsigned long   n_cg;
  unsigned long   i;

  data->ref   = seq_db_new ();

  if (data->verbose)
    g_print (">>> Loading reference %s\n", data->ref_path);
  seq_db_load_fasta (data->ref, data->ref_path, &error);
  if (error)
    {
      g_printerr ("[ERROR] Loading reference failed: %s\n", error->message);
      exit (1);
    }
  n_cg = 0;
  for (i = 0; i < data->ref->total_size; i++)
    if (data->ref->seqs[i] == 'C' || data->ref->seqs[i] == 'G')
      ++n_cg;

  data->meth_index = g_malloc0 (data->ref->total_size * sizeof (*data->meth_index) );
  data->meth_data  = g_malloc0 (n_cg * sizeof (*data->meth_data) );

  n_cg = 0;
  for (i = 0; i < data->ref->total_size; i++)
    if (data->ref->seqs[i] == 'C' || data->ref->seqs[i] == 'G')
      data->meth_index[i] = &data->meth_data[n_cg++];
}

static void
load_previous (CallbackData *data,
               const char   *path)
{
  GIOChannel    *channel;
  SeqDBElement  *elem  = NULL;
  GError        *error = NULL;
  char          *line  = NULL;
  gsize          length;
  gsize          endl;

  if (data->verbose)
    g_print (">>> Loading CG file %s\n", path);

  /* Open */
  channel = g_io_channel_new_file (path, "r", &error);
  if (error)
    {
      g_printerr ("[ERROR] Opening previous results file `%s' failed: %s\n",
                  path, error->message);
      exit (1);
    }

  /* Parse */
  while (G_IO_STATUS_NORMAL == g_io_channel_read_line (channel, &line, &length, &endl, &error))
    {
      line[endl] = '\0';
      if (line[0] == '>')
        {
          elem = g_hash_table_lookup (data->ref->index, line + 1);
          if (!elem)
            g_printerr ("[WARNING] Reference `%s' not found\n", line);
        }
      else if (elem)
        {
          char        **values;
          int           n_fields;
          unsigned long offset;

          values = g_strsplit (line, "\t", -1);

          for (n_fields = 0; values[n_fields]; n_fields++);
          if (n_fields == 3)
            {
              offset = g_ascii_strtoll (values[0], NULL, 10);
              data->meth_index[elem->offset + offset]->n_meth   += g_ascii_strtoll (values[1], NULL, 10);
              data->meth_index[elem->offset + offset]->n_unmeth += g_ascii_strtoll (values[2], NULL, 10);
            }
          else if (n_fields == 4)
            {
              offset = g_ascii_strtoll (values[1], NULL, 10);
              data->meth_index[elem->offset + offset]->n_meth   += g_ascii_strtoll (values[2], NULL, 10);
              data->meth_index[elem->offset + offset]->n_unmeth += g_ascii_strtoll (values[3], NULL, 10);
            }
          else if (n_fields != 1)
            g_printerr ("[WARNING] Could not parse previous run's line\n");

          g_strfreev (values);
        }
      g_free (line);
      line = NULL;
    }
  if (error)
    {
      g_printerr ("[ERROR] Reading previous results file failed: %s\n",
                  error->message);
      exit (1);
    }
  if (line)
    g_free (line);

  /* Close */
  g_io_channel_shutdown (channel, TRUE, &error);
  if (error)
    {
      g_printerr ("[ERROR] Closing previous results file failed: %s\n",
                  error->message);
      g_error_free (error);
    }
  g_io_channel_unref (channel);
}

static void
write_meth (CallbackData *data)
{
  GHashTableIter iter;
  GIOChannel    *channel;
  SeqDBElement  *elem;
  GError        *error      = NULL;
  int            use_stdout = 1;

  if (data->verbose)
    g_print (">>> Writing CG file %s\n", data->output_path);

  /* Open */
  if (data->output_path[0] == '-' && data->output_path[1] == '\0')
    channel = g_io_channel_unix_new (STDOUT_FILENO);
  else
    {
      use_stdout = 0;
      channel = g_io_channel_new_file (data->output_path, "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening output file failed: %s\n", error->message);
          exit (1);
        }
    }

  g_hash_table_iter_init (&iter, data->ref->index);
  while (g_hash_table_iter_next (&iter, NULL, (gpointer*)&elem))
    {
      GString      *buffer;
      MethData    **ref_meth;
      unsigned long i;

      ref_meth = data->meth_index + elem->offset;
      buffer   = g_string_new (NULL);
      g_string_printf (buffer, ">%s\n", elem->name);
      for (i = 0; i < elem->size; i++)
        {
        if (ref_meth[i])
          {
            if (data->print_letter)
              g_string_append_printf (buffer,
                                      "%c\t%ld\t%u\t%u\n",
                                      data->ref->seqs[elem->offset + i],
                                      i,
                                      ref_meth[i]->n_meth,
                                      ref_meth[i]->n_unmeth);
            else
              g_string_append_printf (buffer,
                                      "%ld\t%u\t%u\n",
                                      i,
                                      ref_meth[i]->n_meth,
                                      ref_meth[i]->n_unmeth);
          }
        else if (data->print_all)
          g_string_append_printf (buffer,
                                  "%c\n",
                                  data->ref->seqs[elem->offset + i]);
        }
      g_io_channel_write_chars (channel,
                                buffer->str,
                                buffer->len,
                                NULL,
                                &error);
      if (error)
        {
          g_printerr ("[ERROR] Writing to output file failed: %s\n",
                      error->message);
          g_error_free (error);
          error = NULL;
        }

      g_string_free (buffer, TRUE);
    }

  /* Close */
  if (!use_stdout)
    {
      g_io_channel_shutdown (channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing output file failed: %s\n",
                      error->message);
          g_error_free (error);
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
  if (data->meth_index)
    g_free (data->meth_index);
  if (data->meth_data)
    g_free (data->meth_data);
  seq_db_free (data->ref);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
