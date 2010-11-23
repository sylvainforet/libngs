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

#include "ngs_methylation.h"


typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char              *cg_path;
  char              *ref_path;
  char              *input_path;
  char              *output_path;
  char              *name;

  SeqDB             *ref;
  RefMethCounts     *counts;

  int                from;
  int                to;

  int                verbose;
  int                print_letter;
  int                print_all;
};

static void parse_args     (CallbackData      *data,
                            int               *argc,
                            char            ***argv);

static void load_ref       (CallbackData      *data);

static void cleanup_data   (CallbackData      *data);

static void process_coords (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;

  parse_args (&data, &argc, &argv);

  load_ref (&data);
  process_coords (&data);
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
      {"in",        'i', 0, G_OPTION_ARG_FILENAME, &data->input_path,   "Input file", NULL},
      {"name",      'n', 0, G_OPTION_ARG_STRING,   &data->name,         "Name of the sequence to search", NULL},
      {"from",      'f', 0, G_OPTION_ARG_INT,      &data->from,         "Start position (0-based)", NULL},
      {"to",        't', 0, G_OPTION_ARG_INT,      &data->to,           "End position (exclusive)", NULL},
      {"verbose",   'v', 0, G_OPTION_ARG_NONE,     &data->verbose,      "Verbose output", NULL},
      {"letter",    'l', 0, G_OPTION_ARG_NONE,     &data->print_letter, "Prepend a column with the letter", NULL},
      {"all",       'w', 0, G_OPTION_ARG_NONE,     &data->print_all,    "Prints all positions (implies -l)", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->ref_path     = NULL;
  data->input_path   = NULL;
  data->output_path  = strdup("-");
  data->name         = NULL;
  data->from         = 0;
  data->to           = 0;
  data->verbose      = 0;
  data->print_letter = 0;
  data->print_all    = 0;

  context = g_option_context_new ("FILE - Extracts coordinates from a CG file");
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
  if ((!data->input_path && !data->name) ||
      (data->input_path  &&  data->name))
    {
      g_printerr ("[ERROR] You must either specify an input path with -i or\n"
                  "a name, start and stop with -n, -f and -t respectively\n");
      exit (1);
    }
  if (*argc != 2)
    {
      g_printerr ("[ERROR] You must provide one CG file as argument\n");
      exit (1);
    }
  data->cg_path = (*argv)[1];

  if (data->print_all)
    data->print_letter = 1;
}

static void
load_ref (CallbackData *data)
{
  GError         *error = NULL;

  data->ref   = seq_db_new ();

  if (data->verbose)
    g_print (">>> Loading reference %s\n", data->ref_path);
  seq_db_load_fasta (data->ref,
                     data->ref_path,
                     &error);
  if (error)
    {
      g_printerr ("[ERROR] Loading reference failed: %s\n", error->message);
      exit (1);
    }
  if (data->verbose)
    g_print (">>> Loading CG file: %s\n", data->cg_path);
  data->counts = ref_meth_counts_load (data->ref,
                                       data->cg_path,
                                       &error);
  if (error)
    {
      g_printerr ("[ERROR] failed to load CG file `%s': %s\n",
                  data->cg_path,
                  error->message);
      exit (1);
    }
}


static void
process_coords (CallbackData *data)
{
  GIOChannel *output_channel;
  GError     *error      = NULL;
  int         use_stdout = 1;

  /* Open */
  if (!data->output_path || !*data->output_path || (data->output_path[0] == '-' && data->output_path[1] == '\0'))
    output_channel = g_io_channel_unix_new (STDOUT_FILENO);
  else
    {
      use_stdout     = 0;
      output_channel = g_io_channel_new_file (data->output_path, "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] failed to open output file `%s': %s\n",
                      data->output_path,
                      error->message);
          exit (1);
        }
    }

  if (data->name)
    {
      ref_meth_counts_write_segment (data->counts,
                                     data->ref,
                                     output_channel,
                                     data->name,
                                     MAX (data->from, 0),
                                     MAX (data->to, 0) ,
                                     data->print_letter,
                                     data->print_all,
                                     &error);
      if (error)
        {
          g_printerr ("[ERROR] Problem while writing segment: %s\n",
                      error->message);
          g_error_free (error);
        }
    }
  else
    {
      GIOChannel  *input_channel;
      char       **fields;
      char        *line;
      gsize        length;
      gsize        endl;


      /* Open */
      input_channel = g_io_channel_new_file (data->input_path, "r", &error);
      if (error)
        {
          g_printerr ("[ERROR] failed to open input file `%s': %s\n",
                      data->input_path,
                      error->message);
          exit (1);
        }

      while (G_IO_STATUS_NORMAL == g_io_channel_read_line (input_channel, &line, &length, &endl, &error))
        {
          int n_fields;

          line[endl] = '\0';
          fields     = g_strsplit_set (line, " \t", 0);
          n_fields   = -1;
          while (fields[++n_fields]);
          if (n_fields < 3)
            g_printerr ("[WARNING] Could not parse coordinates line: %s\n", line);
          else
            {
              long from;
              long to;

              from = g_ascii_strtoll (fields[1], NULL, 10);
              to   = g_ascii_strtoll (fields[2], NULL, 10);
              ref_meth_counts_write_segment (data->counts,
                                             data->ref,
                                             output_channel,
                                             fields[0],
                                             MAX (from, 0),
                                             MAX (to, 0) ,
                                             data->print_letter,
                                             data->print_all,
                                             &error);
            }
          g_strfreev (fields);
          g_free (line);
        }

      /* Close */
      g_io_channel_shutdown (input_channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing input file `%s' failed: %s\n",
                      data->input_path,
                      error->message);
          g_error_free (error);
        }
      g_io_channel_unref (input_channel);
    }

  /* Close */
  if (!use_stdout)
    {
      g_io_channel_shutdown (output_channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing output file `%s' failed: %s\n",
                      data->output_path,
                      error->message);
          g_error_free (error);
        }
    }
  g_io_channel_unref (output_channel);
}

static void
cleanup_data (CallbackData *data)
{
  if (!data)
    return;
  if (data->ref_path)
    g_free (data->ref_path);
  if (data->input_path)
    g_free (data->input_path);
  if (data->output_path)
    g_free (data->output_path);
  if (data->name)
    g_free (data->name);
  if (data->counts)
    ref_meth_counts_destroy (data->counts);
  seq_db_free (data->ref);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
