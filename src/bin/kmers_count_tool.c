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
#include <unistd.h>
#include <string.h>

#include "ngs_kmerhash.h"

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  KmerHashTable     *hash_table;

  char             **in_txt;
  char             **in_bin;
  char              *out;

  unsigned int       k;

  int                verbose;
  int                binary;
};

static void parse_args    (CallbackData      *data,
                           int               *argc,
                           char            ***argv);

static void cleanup       (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;
  char        **tmp;
  GError       *error = NULL;

  parse_args (&data, &argc, &argv);

  if (data.in_txt)
    for (tmp = data.in_txt; *tmp != NULL; tmp++)
      {
        if (data.verbose)
          g_printerr ("Loading kmer counts from text file: %s\n", *tmp);
        kmer_hash_table_load (data.hash_table, *tmp, 0, &error);
        if (error)
          {
            g_printerr ("[ERROR] Problem while loading %s: %s\n",
                        *tmp, error->message);
            exit (1);
          }
      }
  if (data.in_bin)
    for (tmp = data.in_bin; *tmp != NULL; tmp++)
      {
        if (data.verbose)
          g_printerr ("Loading kmer counts from binary file: %s\n", *tmp);
        kmer_hash_table_load (data.hash_table, *tmp, 1, &error);
        if (error)
          {
            g_printerr ("[ERROR] Problem while loading %s: %s\n",
                        *tmp, error->message);
            exit (1);
          }
      }
  if (data.verbose)
    g_printerr ("Writing results to %s\n", data.out);
  kmer_hash_table_print (data.hash_table,
                         data.out,
                         0,
                         &error);
  if (error)
    {
      g_printerr ("[ERROR] Problem while writing %s: %s\n",
                  data.out, error->message);
      exit (1);
    }

  cleanup (&data);

  return 0;
}

static void
parse_args (CallbackData      *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      {"kmer",     'k', 0, G_OPTION_ARG_INT,            &data->k,        "K-mer size",                   NULL},
      {"verbose",  'v', 0, G_OPTION_ARG_NONE,           &data->verbose,  "Verbose output",               NULL},
      {"binary",   'u', 0, G_OPTION_ARG_NONE,           &data->binary,   "Output in binary format",      NULL},
      {"out",      'o', 0, G_OPTION_ARG_FILENAME,       &data->out,      "Output file",                  NULL},
      {"intxt",    'a', 0, G_OPTION_ARG_FILENAME_ARRAY, &data->in_txt,   "Input text kmer count file",   NULL},
      {"inbin",    'b', 0, G_OPTION_ARG_FILENAME_ARRAY, &data->in_bin,   "Input binary kmer count file", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->k        = 0;
  data->verbose  = 0;
  data->binary   = 0;
  data->out      = NULL;
  data->in_txt   = NULL;
  data->in_bin   = NULL;

  context = g_option_context_new ("- Utility to manipulate kmers count files");
  g_option_context_add_main_entries (context, entries, NULL);
  if (!g_option_context_parse (context, argc, argv, &error))
    {
      g_printerr ("[ERROR] Option parsing failed: %s\n", error->message);
      exit (1);
    }
  g_option_context_free (context);

  if (data->k < 1)
    {
      g_printerr ("[ERROR] The kmer size must be larger than 0\n");
      exit (1);
    }
  if (!data->out)
    data->out = g_strdup ("-");
  if ((data->in_txt == NULL || *data->in_txt == NULL) &&
      (data->in_bin == NULL || *data->in_bin == NULL))
    {
      g_printerr ("No input file provided\n");
      exit (1);
    }
  data->hash_table = kmer_hash_table_new (data->k);
}

static void
cleanup (CallbackData *data)
{
  if (data->in_txt)
    g_strfreev (data->in_txt);
  if (data->in_bin)
    g_strfreev (data->in_bin);
  if (data->out)
    g_free (data->out);
  if (data->hash_table)
    kmer_hash_table_destroy (data->hash_table);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
