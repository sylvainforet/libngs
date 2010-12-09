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

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char       *input_path;
  char       *output_path;
  GIOChannel *output_channel;

  int         start;
  int         end;
  int         tot_trim;
  int         qual;
  int         len;
  int         win;
  int         qwin;
  int         keep;
  int         sanger;
  int         qual0;

  int         use_stdout: 1;
};

static void parse_args  (CallbackData      *data,
                         int               *argc,
                         char            ***argv);

static int  iter_func   (FastqSeq          *fastq,
                         CallbackData      *data);

static int  write_fastq (GIOChannel        *channel,
                         char              *name,
                         char              *seq,
                         char              *qual,
                         int                start,
                         int                end);

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
      g_printerr ("[ERROR] Iterating sequences failed: %s\n",
                  error->message);
      g_error_free (error);
      error = NULL;
      return 1;
    }
  if (data.output_channel)
    {
      if (!data.use_stdout)
        {
          g_io_channel_shutdown (data.output_channel, TRUE, &error);
          if (error)
            {
              g_printerr ("[ERROR] Closing output file failed: %s\n",
                          error->message);
              g_error_free (error);
              error = NULL;
            }
        }
      g_io_channel_unref (data.output_channel);
    }

  return 0;
}

static void
parse_args (CallbackData      *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      {"out",    'o', 0, G_OPTION_ARG_FILENAME, &data->output_path, "Output file", NULL},
      {"start",  's', 0, G_OPTION_ARG_INT,  &data->start,  "Number of positions to trim at the start", NULL},
      {"end",    'e', 0, G_OPTION_ARG_INT,  &data->end,    "Number of positions to trim at the end", NULL},
      {"qual",   'q', 0, G_OPTION_ARG_INT,  &data->qual,   "Minimum quality", NULL},
      {"len",    'l', 0, G_OPTION_ARG_INT,  &data->len,    "Minimum length", NULL},
      {"win",    'w', 0, G_OPTION_ARG_INT,  &data->win,    "Size of the slidding window", NULL},
      {"qwin",   'n', 0, G_OPTION_ARG_INT,  &data->qwin,   "Minimum sliding window mean quality", NULL},
      {"keep",   'k', 0, G_OPTION_ARG_NONE, &data->keep,   "Keep an pseudo-entry for too small reads", NULL},
      {"sanger", 'g', 0, G_OPTION_ARG_NONE, &data->sanger, "Qualities are in sanger format", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->output_path    = "-";
  data->output_channel = NULL;
  data->use_stdout     = 1;
  data->start          = 0;
  data->end            = 0;
  data->qual           = 0;
  data->len            = 0;
  data->win            = 0;
  data->qwin           = 0;
  data->keep           = 0;
  data->sanger         = 0;
  data->qual0          = 0;

  context = g_option_context_new ("FILE - trims fastq reads");
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
  if (data->len == 0 && data->win > 0)
    data->len = data->win;
  if (data->win > data->len)
    {
      g_printerr ("[ERROR] Size of sliding window (%d) exceeds minimum length (%d)\n",
                  data->win, data->len);
      exit (1);
    }
  if (data->qwin > 0 && data->qual > data->qwin)
    {
      g_printerr ("[ERROR] Minimum quality (%d) exceeds sliding window quality (%d)\n",
                  data->qual, data->qwin);
      exit (1);
    }
  data->tot_trim = data->start + data->end;

  data->qual0 = FASTQ_QUAL_0;
  if (data->sanger)
    data->qual0 = FASTQ_QUAL_0_SANGER;

  if (data->output_path[0] == '-' && data->output_path[1] == '\0')
    data->output_channel = g_io_channel_unix_new (STDOUT_FILENO);
  else
    {
      data->use_stdout     = 0;
      data->output_channel = g_io_channel_new_file (data->output_path, "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening output file failed: %s\n", error->message);
          exit (1);
        }
    }
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  int start;
  int end;

  if (data->tot_trim > fastq->size)
    {
      g_printerr ("[ERROR] trimming more than sequence length "
                  "(trim: %d - length: %d)\n",
                  data->tot_trim, fastq->size);
      exit (1);
    }

  end     = fastq->size - data->end;
  start   = data->start;
  /* Low quality bases */
  if (data->qual > 0)
    {
      int max_idx = start;
      int max_len = 0;
      int tmp_idx = start;
      int i;
      for (i = start; i < end; i++)
        if (fastq->qual[i] - data->qual0 < data->qual)
          {
            if (i - tmp_idx > max_len)
              {
                max_idx = tmp_idx;
                max_len = i - tmp_idx;
              }
            tmp_idx = i + 1;
          }
      if (i - tmp_idx > max_len)
        {
          max_idx = tmp_idx;
          max_len = i - tmp_idx;
        }
      start = max_idx;
      end   = max_idx + max_len;
    }
  /* TODO
   * This code for sliding window is not optimised, this could use some work.
   * I should also make sure that everything is consistent regarding the
   * various combinations of options.
   */
  /* Sliding window */
  if (data->qwin > 0 && data->win > 0 && end - start >= data->len)
    {
      int tmp_len  = 0;
      int max_len  = 0;
      int win_qual = 0;
      int good_win = 0;
      int max_idx  = start;
      int i;
      /* Initialise window */
      for (i = start + data->win - 1; i >= start; i--)
        win_qual += fastq->qual[i] - data->qual0;
      for (i = start; i <= end - data->win; i++)
        {
          if (i > start)
            {
              win_qual -= fastq->qual[i - 1] - data->qual0;
              win_qual += fastq->qual[i + data->win - 1] - data->qual0;
            }
          if ((win_qual / data->win) < data->qwin) /* Bad quality window */
            {
              if (good_win)
                {
                  good_win = 0;
                  if (tmp_len > max_len)
                    {
                      max_len = tmp_len;
                      max_idx = i + data->win - max_len - 1;
                    }
                  tmp_len = 0;
                }
            }
          else /* Good quality window */
            {
              if (good_win)
                tmp_len++;
              else
                {
                  tmp_len  = data->win;
                  good_win = 1;
                }
            }
        }
      if (good_win && tmp_len > max_len)
        {
          max_len = tmp_len;
          max_idx = i + data->win - max_len - 1;
        }
      start = max_idx;
      end   = max_idx + max_len;
    }

  if (end <= start || end - start < data->len)
    {
      if (data->keep)
        return write_fastq (data->output_channel,
                            fastq->name,
                            "N", "2",  0, 1);
      return 1;
    }
  return write_fastq (data->output_channel,
                      fastq->name,
                      fastq->seq,
                      fastq->qual,
                      start,
                      end);
}

static int
write_fastq (GIOChannel *channel,
             char       *name,
             char       *seq,
             char       *qual,
             int         start,
             int         end)
{
  GError  *error = NULL;
  GString *buffer;
  int      ret = 1;

  buffer = g_string_sized_new (512);

  g_string_append_printf (buffer, "@%s\n", name);

  buffer = g_string_append_len (buffer, seq + start, end - start);
  buffer = g_string_append_c (buffer, '\n');

  g_string_append_printf (buffer, "+%s\n", name);

  buffer = g_string_append_len (buffer, qual + start, end - start);
  buffer = g_string_append_c (buffer, '\n');

  g_io_channel_write_chars (channel,
                            buffer->str,
                            -1,
                            NULL,
                            &error);
  if (error)
    {
      g_printerr ("[ERROR] Writing trimmed sequence failed: %s\n",
                  error->message);
      g_error_free (error);
      error = NULL;
      ret   = 0;
    }
  g_string_free (buffer, TRUE);

  return ret;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
