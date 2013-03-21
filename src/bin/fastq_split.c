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

#include "ngs_fastq.h"

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char       *input_path;
  GIOChannel *output_channel;
  GString    *buffer;

  char       *prefix;
  int         chunks;
  int         count;
  int         reads;
};

static int  iter_func  (FastqSeq       *fastq,
                        CallbackData   *data);

static void parse_args (CallbackData   *data,
                        int            *argc,
                        char         ***argv);

int
main (int    argc,
      char **argv)
{
  CallbackData data;
  GError      *error = NULL;

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
    }

  /* Cleanup */
  if (data.output_channel)
    {
      g_io_channel_shutdown (data.output_channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing sequence output file failed: %s\n", error->message);
          g_error_free (error);
          error = NULL;
        }
      g_io_channel_unref (data.output_channel);
    }
  g_string_free (data.buffer, TRUE);

  return 0;
}

static void
parse_args (CallbackData   *data,
            int            *argc,
            char         ***argv)
{
  GOptionEntry entries[] =
    {
      {"prefix", 'p', 0, G_OPTION_ARG_STRING, &data->prefix, "Prefix for the output chunks", NULL},
      {"reads",  'r', 0, G_OPTION_ARG_INT,    &data->reads,  "Number of reads per chunk",    NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->output_channel = NULL;
  data->prefix         = NULL;
  data->chunks         = 0;
  data->count          = 0;
  data->reads          = G_MAXINT;

  context = g_option_context_new ("FILE - Converts between various flavours of fastq formats");
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

  if (!data->prefix)
    data->prefix = "chunk";

  data->count = data->reads + 1;

  data->buffer = g_string_sized_new (1024);
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  GError  *error = NULL;
  int      ret   = 1;

  if (data->count > data->reads)
    {
      char *path;

      /* Close current output */
      if (data->output_channel)
        {
          g_io_channel_shutdown (data->output_channel, TRUE, &error);
          if (error)
            {
              g_printerr ("[ERROR] Closing sequence output file failed: %s\n", error->message);
              g_error_free (error);
              error = NULL;
            }
          g_io_channel_unref (data->output_channel);
        }

      /* Open new output */
      path = g_strdup_printf ("%s_%06d", data->prefix, data->chunks);
      data->output_channel = g_io_channel_new_file (path, "w", &error);
      g_free (path);
      if (error)
        {
          data->output_channel = NULL;
          g_printerr ("[ERROR] Opening sequence output file failed: %s\n", error->message);
          g_error_free (error);
          return 0;
        }

      data->chunks += 1;
      data->count = 1;
    }

  /* Write current read */
  fastq_write (data->output_channel,
               data->buffer,
               fastq->name,
               fastq->seq,
               fastq->qual,
               &error);
  if (error)
    {
      g_printerr ("[ERROR] Writing sequence failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
      ret   = 0;
    }

  data->count += 1;

  return ret;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
*/
