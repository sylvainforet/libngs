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
  unsigned long int  qual[128];

  char              *input_path;

  unsigned int       all;
};

static void parse_args    (CallbackData      *data,
                           int               *argc,
                           char            ***argv);

static int  iter_func_all (FastqSeq          *fastq,
                           CallbackData      *data);

static int  iter_func     (FastqSeq          *fastq,
                           CallbackData      *data);

static void init_qual     (CallbackData      *data);

static void print_qual    (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;
  GError       *error = NULL;

  parse_args (&data, &argc, &argv);

  if (data.all)
    iter_fastq (data.input_path,
                (FastqIterFunc)iter_func_all,
                &data,
                &error);
  else
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

  if (!data.all)
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
      {"all",    'a', 0, G_OPTION_ARG_NONE, &data->all,    "Do not bin the qualities, print quality for each read", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->all = 0;

  context = g_option_context_new ("FILE - distribution of reads qualities, or lists all reads qualities (-a option)");
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

  if (!data->all)
    init_qual (data);
}

static int
iter_func_all (FastqSeq     *fastq,
               CallbackData *data)
{
  int i;
  long read_qual = 0;

  for (i = 0; i < fastq->size; i++)
    read_qual += fastq->qual[i] - fastq_qual0;

  g_print ("%f\n", ((float)read_qual) / fastq->size);

  return 1;
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  int i;
  long read_qual = 0;
  long read_mean = 0;

  for (i = 0; i < fastq->size; i++)
    read_qual += fastq->qual[i] - fastq_qual0;

  read_mean = 0.5 + ((float)read_qual) / fastq->size;
  data->qual[read_mean]++;

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
