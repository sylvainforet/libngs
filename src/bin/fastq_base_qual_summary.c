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

#include "ngs_fastq.h"


#define N_QUAL 60

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  unsigned long int  qual[N_QUAL];

  char              *input_path;
};

static void parse_args (CallbackData      *data,
                        int               *argc,
                        char            ***argv);

static int  iter_func  (FastqSeq          *fastq,
                        CallbackData      *data);

static void init_qual  (CallbackData      *data);

static void print_qual (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;
  GError       *error = NULL;

  parse_args (&data, &argc, &argv);

  iter_fastq (data.input_path,
              (FastqIterFunc)iter_func,
              &data,
              &error);

  if (error)
    {
      g_printerr ("[ERROR] Iterating sequences failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
      return 1;
    }

  print_qual (&data);

  return 0;
}

static void
parse_args (CallbackData      *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      /*{"fast", 'f', 0, G_OPTION_ARG_NONE, &data->fast, "Faster and less robust method (use at your own risk)", NULL},*/
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  context = g_option_context_new ("FILE - Computes distribution of the base qualities");
  g_option_context_add_group (context, get_fastq_option_group ());
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

  init_qual (data);
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  int i;

  for (i = 0; i < fastq->size; i++)
    data->qual[fastq->qual[i] - fastq_qual0]++;

  return 1;
}

static void
init_qual (CallbackData *data)
{
  int i;

  for (i = 0; i < N_QUAL; i++)
    data->qual[i] = 0;
}

static void
print_qual (CallbackData *data)
{
  int i;

  for (i = 0; i < N_QUAL; i++)
    g_print ("%d\t%ld\n",
             i, data->qual[i]);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
