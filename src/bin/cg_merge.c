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
  char              *ref_path;
  char              *output_path;

  SeqDB             *ref;
  RefMethCounts     *counts;

  int                verbose;
  int                print_letter;
  int                print_all;
};

static void parse_args    (CallbackData      *data,
                           int               *argc,
                           char            ***argv);

static void load_ref      (CallbackData      *data);

static void cleanup_data  (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;
  GError       *error = NULL;
  int           i;

  parse_args (&data, &argc, &argv);

  load_ref (&data);
  for (i = 1; i < argc; i++)
    {
      if (data.verbose)
        g_print (">>> Loading meth file: %s\n", argv[i]);
      ref_meth_counts_add_path (data.counts,
                                data.ref,
                                argv[i],
                                &error);
      if (error)
        {
          g_printerr ("[ERROR] failed to load meth file `%s': %s\n",
                      argv[i],
                      error->message);
          exit (1);
        }
    }
  if (data.verbose)
    g_print (">>> Writing meth file: %s\n", data.output_path);
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
      {"reference", 'r', 0, G_OPTION_ARG_FILENAME, &data->ref_path,     "Reference genome file", NULL},
      {"out",       'o', 0, G_OPTION_ARG_FILENAME, &data->output_path,  "Output file", NULL},
      {"verbose",   'v', 0, G_OPTION_ARG_NONE,     &data->verbose,      "Verbose output", NULL},
      {"letter",    'l', 0, G_OPTION_ARG_NONE,     &data->print_letter, "Prepend a column with the letter", NULL},
      {"all",       'w', 0, G_OPTION_ARG_NONE,     &data->print_all,    "Prints all positions (implies -l)", NULL},
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
  data->counts = ref_meth_counts_create (data->ref);
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
