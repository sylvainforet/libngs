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

#include "ngs_fasta.h"

typedef struct  _TestFastaData TestFastaData;

struct _TestFastaData
{
  char *input_path;
  int   width;
};

static void parse_args (TestFastaData    *data,
                        int               *argc,
                        char            ***argv);

int
main (int    argc,
      char **argv)
{
  TestFastaData  data;
  FastaIter     *iter;
  FastaSeq      *seq;
  GError        *error = NULL;

  parse_args (&data, &argc, &argv);
  iter = fasta_iter_new (data.input_path, &error);
  if (error)
    {
      g_printerr ("[ERROR] Failed to load input file `%s': %s\n",
                  data.input_path,
                  error->message);
      exit (1);
    }
  for (seq = fasta_iter_next (iter) ; seq != NULL; seq = fasta_iter_next (iter))
    {
      char *tmp;

      tmp = fasta_write_to_buffer (seq, 50);
      g_print ("%s\n", tmp);
      g_free (tmp);
    }

  fasta_iter_free (iter);

  return 0;
}

static void
parse_args (TestFastaData    *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      /* {"width", 'w', 0, G_OPTION_ARG_INT, &data->width, "Line width", NULL}, */
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->width = 50;

  context = g_option_context_new ("FILE - Reads a fasta file and spits it back out");
  g_option_context_add_group (context, get_fasta_option_group ());
  g_option_context_add_main_entries (context, entries, NULL);
  if (!g_option_context_parse (context, argc, argv, &error))
    {
      g_printerr ("[ERROR] Option parsing failed: %s\n", error->message);
      exit (1);
    }
  g_option_context_free (context);

  if (*argc < 2)
    {
      g_printerr ("[ERROR] No input file provided\n");
      exit (1);
    }
  data->input_path = (*argv)[1];
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
