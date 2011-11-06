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
  char       *output_path;
  GIOChannel *output_channel;

  char       *qual0;
  int         qual_max;
  int         delta_qual;
  int         qual_name;
  int         old_pairs;

  int         use_stdout;

  char        qual0c;
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

  if (data.output_channel)
    {
      if (!data.use_stdout)
        {
          g_io_channel_shutdown (data.output_channel, TRUE, &error);
          if (error)
            {
              g_printerr ("[ERROR] Closing sequence output file failed: %s\n", error->message);
              g_error_free (error);
              error = NULL;
            }
        }
      g_io_channel_unref (data.output_channel);
    }

  return 0;
}

static void
parse_args (CallbackData   *data,
            int            *argc,
            char         ***argv)
{
  GOptionEntry entries[] =
    {
      {"out"  ,        'o', 0, G_OPTION_ARG_FILENAME, &data->output_path, "Output file",                                 NULL},
      {"qual0",        'q', 0, G_OPTION_ARG_STRING,   &data->qual0,       "The character for 0 quality value",           NULL},
      {"qual_max",     'Q', 0, G_OPTION_ARG_INT,      &data->qual_max,    "Maximum quality value",                       NULL},
      {"qual_name",    'n', 0, G_OPTION_ARG_NONE,     &data->qual_name,   "Write the full name before the quality line", NULL},
      {"old_pairs",    'p', 0, G_OPTION_ARG_NONE,     &data->old_pairs,   "Convert new pair format to old format",       NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->output_path    = NULL;
  data->output_channel = NULL;
  data->qual0          = NULL;
  data->qual0c         = fastq_qual0;
  data->qual_max       = -1;
  data->qual_name      = 0;
  data->old_pairs      = 0;
  data->use_stdout     = 0;

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

  if (data->qual0 != NULL)
    {
      if (strlen (data->qual0) != 1)
        {
          g_printerr ("[ERROR] The `-qual0/-q' option expects a single character\n");
          exit (1);
        }
      /* TODO some quality checks on the new value */
      data->qual0c = data->qual0[0];
    }
  data->delta_qual = data->qual0c - fastq_qual0;

  if (!data->output_path)
    data->output_path = g_strdup ("-");

  if (data->output_path[0] == '-' && data->output_path[1] == '\0')
    {
      data->use_stdout     = 1;
      data->output_channel = g_io_channel_unix_new (STDOUT_FILENO);
    }
  else
    {
      data->output_channel = g_io_channel_new_file (data->output_path, "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening sequence output file failed: %s\n", error->message);
          exit (1);
        }
    }
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  GString *buffer;
  GError  *error = NULL;
  int      ret   = 1;

  /* A rather hugly hack to change the name in place */
  if (data->old_pairs)
    {
      int       i;
      int       n_separators = 0;
      int       pair_idx     = -1;
      const int l            = strlen(fastq->name);

      for (i = 0; i < l; i++)
        {
          if (fastq->name[i] == ':')
            ++n_separators;
          if (n_separators == 7)
            {
              if (fastq->name[i - 1] == '1')
                pair_idx = '1';
              else if (fastq->name[i - 1] == '2')
                pair_idx = '2';
              break;
            }
        }
      if (pair_idx > 0)
        {
          for (; i < l; i++)
            fastq->name[i - 2] = fastq->name[i];

          fastq->name[l - 2] = '/';
          fastq->name[l - 1] = pair_idx;
        }
    }
  
  /* Quality convertion */
  if (data->delta_qual != 0)
    {
      int i;

      for (i = 0; i < fastq->size; i++)
        fastq->qual[i] += data->delta_qual;
    }
  /* Maximum quality, to keep allpaths happy */
  if (data->qual_max >= 0)
    {
      int i;

      for (i = 0; i < fastq->size; i++)
        if (fastq->qual[i] > data->qual_max)
          fastq->qual[i] = data->qual_max;
    }

  buffer = g_string_sized_new (512);
  buffer = g_string_append_c (buffer, '@');
  buffer = g_string_append (buffer, fastq->name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append (buffer, fastq->seq);
  buffer = g_string_append (buffer, "\n+");
  if (data->qual_name)
    buffer = g_string_append (buffer, fastq->name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append (buffer, fastq->qual);
  buffer = g_string_append_c (buffer, '\n');

  g_io_channel_write_chars (data->output_channel,
                            buffer->str,
                            -1,
                            NULL,
                            &error);
  if (error)
    {
      g_printerr ("[ERROR] Writing sequence failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
      ret   = 0;
    }
  g_string_free (buffer, TRUE);

  return ret;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
*/
