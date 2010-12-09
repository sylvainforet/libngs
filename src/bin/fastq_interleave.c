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

#include "ngs_fastq.h"
#include "ngs_utils.h"

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char       *input_path1;
  char       *input_path2;
  char       *output_path;
  char       *output_single;
  GIOChannel *output_channel;
  GIOChannel *single_channel;

  int         min_size;
  int         append;
  int         append_single;

  int         use_stdout: 1;
};

static void parse_args  (CallbackData   *data,
                         int            *argc,
                         char         ***argv);

static void write_fastq (FastqSeq       *fastq,
                         GString        *buffer,
                         GIOChannel     *channel);

int
main (int    argc,
      char **argv)
{
  CallbackData data;
  FastqSeq    *seq1;
  FastqSeq    *seq2;
  FastqIter   *iter1;
  FastqIter   *iter2;
  GString     *buffer;
  GError      *error = NULL;

  parse_args (&data, &argc, &argv);

  iter1 = fastq_iter_new (data.input_path1, &error);
  if (error)
    {
      g_printerr ("[ERROR] Opening %s failed: %s\n", data.input_path1, error->message);
      exit (1);
    }
  iter2 = fastq_iter_new (data.input_path2, &error);
  if (error)
    {
      g_printerr ("[ERROR] Opening %s failed: %s\n", data.input_path2, error->message);
      exit (1);
    }


  buffer = g_string_sized_new (1024);
  seq1 = fastq_iter_next (iter1);
  seq2 = fastq_iter_next (iter2);
  while (seq1 != NULL && seq2 != NULL)
    {
      if (data.min_size > 0)
        {
          if (seq1->size >= data.min_size && seq2->size >= data.min_size)
            {
              write_fastq (seq1, buffer, data.output_channel);
              write_fastq (seq2, buffer, data.output_channel);
            }
          else
            {
              if (seq1->size >= data.min_size && data.single_channel)
                write_fastq (seq1, buffer, data.single_channel);
              if (seq2->size >= data.min_size && data.single_channel)
                write_fastq (seq2, buffer, data.single_channel);
            }
        }
      else
        {
          write_fastq (seq1, buffer, data.output_channel);
          write_fastq (seq2, buffer, data.output_channel);
        }
      seq1 = fastq_iter_next (iter1);
      seq2 = fastq_iter_next (iter2);
    }
  fastq_iter_free (iter1);
  fastq_iter_free (iter2);
  g_string_free (buffer, TRUE);

  if (data.output_channel)
    {
      if (!data.use_stdout)
        {
          g_io_channel_shutdown (data.output_channel, TRUE, &error);
          if (error)
            {
              g_printerr ("[ERROR] Closing output file failed: %s\n", error->message);
              g_error_free (error);
              error = NULL;
            }
        }
      g_io_channel_unref (data.output_channel);
    }
  if (data.output_path)
    g_free (data.output_path);
  if (data.single_channel)
    {
      g_io_channel_shutdown (data.single_channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing single reads output file failed: %s\n", error->message);
          g_error_free (error);
          error = NULL;
        }
      g_io_channel_unref (data.single_channel);
    }
  if (data.output_single)
    g_free (data.output_single);

  return 0;
}

static void
parse_args (CallbackData   *data,
            int            *argc,
            char         ***argv)
{
  GOptionEntry entries[] =
    {
      {"output",       'o', 0, G_OPTION_ARG_FILENAME, &data->output_path,   "Output path", NULL},
      {"single",       's', 0, G_OPTION_ARG_FILENAME, &data->output_single, "File for single reads", NULL},
      {"minsize",      'm', 0, G_OPTION_ARG_INT,      &data->min_size,      "Min read size", NULL},
      {"append",       'a', 0, G_OPTION_ARG_NONE,     &data->append,        "Append to output", NULL},
      {"appendsingle", 'b', 0, G_OPTION_ARG_NONE,     &data->append_single, "Append to single reads output", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->output_path    = NULL;
  data->output_single  = NULL;
  data->output_channel = NULL;
  data->single_channel = NULL;
  data->append         = 0;
  data->append_single  = 0;
  data->use_stdout     = 1;
  data->min_size       = 0;

  context = g_option_context_new ("FILE1 FILE2 - interleaves the sequences from two fastq files");
  g_option_context_add_group (context, get_fastq_option_group ());
  g_option_context_add_main_entries (context, entries, NULL);
  if (!g_option_context_parse (context, argc, argv, &error))
    {
      g_printerr ("[ERROR] Option parsing failed: %s\n", error->message);
      exit (1);
    }
  g_option_context_free (context);

  if (*argc < 3)
    {
      g_printerr ("[ERROR] No input file provided\n");
      exit (1);
    }
  data->input_path1 = (*argv)[1];
  data->input_path2 = (*argv)[2];

  if (!data->output_path)
    data->output_path = g_strdup ("-");
  if (data->output_path[0] == '-' && data->output_path[1] == '\0')
    {
      data->use_stdout      = 1;
      data->output_channel = g_io_channel_unix_new (STDOUT_FILENO);
    }
  else
    {
      data->use_stdout      = 0;
      data->output_channel = g_io_channel_new_file (data->output_path,
                                                    data->append ? "a" : "w",
                                                    &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening sequence output file failed: %s\n", error->message);
          exit (1);
        }
    }
  g_io_channel_set_encoding (data->output_channel, NULL, NULL);
  g_io_channel_set_buffered (data->output_channel, TRUE);
  g_io_channel_set_buffer_size (data->output_channel, 4096 * 128);
  if (data->output_single)
    {
      data->single_channel = g_io_channel_new_file (data->output_single,
                                                    data->append_single ? "a" : "w",
                                                    &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening single reads output file failed: %s\n", error->message);
          exit (1);
        }
    }
  g_io_channel_set_encoding (data->single_channel, NULL, NULL);
  g_io_channel_set_buffered (data->single_channel, TRUE);
  g_io_channel_set_buffer_size (data->single_channel, 4096 * 128);
}

static void
write_fastq (FastqSeq   *fastq,
             GString    *buffer,
             GIOChannel *channel)
{
  GError  *error = NULL;

  g_string_truncate (buffer, 0);

  /* Sequence */
  buffer = g_string_append_c (buffer, '@');
  buffer = g_string_append (buffer, fastq->name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append (buffer, fastq->seq);
  buffer = g_string_append_c (buffer, '\n');

  /* Quality */
  buffer = g_string_append_c (buffer, '+');
  buffer = g_string_append (buffer, fastq->name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append (buffer, fastq->qual);
  buffer = g_string_append_c (buffer, '\n');

  g_io_channel_write_chars (channel,
                            buffer->str,
                            -1,
                            NULL,
                            &error);
  if (error)
    {
      g_printerr ("[ERROR] Writing sequence failed: %s\n", error->message);
      exit (1);
    }
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
