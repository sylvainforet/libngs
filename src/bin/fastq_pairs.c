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
  char       *output_path1;
  char       *output_path2;
  char       *output_single;
  GIOChannel *output_channel1;
  GIOChannel *output_channel2;
  GIOChannel *single_channel;

  int         min_size;
  int         append;
  int         append_single;

  int         one_input;
  int         one_output;
  int         use_stdout1;
  int         use_stdout2;
  int         use_stdout_single;
};

static void parse_args   (CallbackData   *data,
                          int            *argc,
                          char         ***argv);

static void cleanup      (CallbackData   *data);

static void open_output  (GIOChannel **channel,
                          int         *use_stdout,
                          const char  *path,
                          const char  *mode,
                          const char  *error_message);

static void close_output (GIOChannel **channel,
                          int          use_stdout,
                          const char * error_message);

int
main (int    argc,
      char **argv)
{
  CallbackData data;
  FastqSeq    *seq1;
  FastqSeq    *seq2;
  FastqSeq    *tmp;
  FastqIter   *iter1;
  FastqIter   *iter2;
  GString     *buffer;
  GError      *error = NULL;

  parse_args (&data, &argc, &argv);

  iter1 = fastq_iter_new (data.input_path1, &error);
  if (error)
    {
      g_printerr ("[ERROR] Opening %s failed: %s\n",
                  data.input_path1, error->message);
      exit (1);
    }
  if (!data.one_input)
    {
      iter2 = fastq_iter_new (data.input_path2, &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening %s failed: %s\n",
                      data.input_path2, error->message);
          exit (1);
        }
    }
  else
    iter2 = iter1;

  buffer = g_string_sized_new (1024);
  tmp    = fastq_iter_next (iter1);
  if (data.one_input)
    seq1 = fastq_seq_copy (tmp);
  else
    seq1 = tmp;
  seq2   = fastq_iter_next (iter2);
  while (seq1 != NULL && seq2 != NULL)
    {
      if (data.min_size > 0)
        {
          if (seq1->size >= data.min_size && seq2->size >= data.min_size)
            {
              fastq_write (data.output_channel1,
                           buffer,
                           seq1->name,
                           seq1->seq,
                           seq1->qual,
                           &error);
              if (error != NULL)
                break;
              fastq_write (data.output_channel2,
                           buffer,
                           seq2->name,
                           seq2->seq,
                           seq2->qual,
                           &error);
              if (error != NULL)
                break;
            }
          else if (data.single_channel != NULL)
            {
              if (seq1->size >= data.min_size)
                {
                  fastq_write (data.single_channel,
                               buffer,
                               seq1->name,
                               seq1->seq,
                               seq1->qual,
                               &error);
                  if (error != NULL)
                    break;
                }
              if (seq2->size >= data.min_size)
                {
                  fastq_write (data.single_channel,
                               buffer,
                               seq2->name,
                               seq2->seq,
                               seq2->qual,
                               &error);
                  if (error != NULL)
                    break;
                }
            }
        }
      else
        {
          fastq_write (data.output_channel1,
                       buffer,
                       seq1->name,
                       seq1->seq,
                       seq1->qual,
                       &error);
          if (error != NULL)
            break;
          fastq_write (data.output_channel2,
                       buffer,
                       seq2->name,
                       seq2->seq,
                       seq2->qual,
                       &error);
          if (error != NULL)
            break;
        }
      if (data.one_input)
        fastq_seq_free (seq1);
      tmp = fastq_iter_next (iter1);
      if (data.one_input)
        seq1 = fastq_seq_copy (tmp);
      else
        seq1 = tmp;
      seq2 = fastq_iter_next (iter2);
    }
  /* In case seq1 was allocated, but the loop was never called */
  if (data.one_input)
    fastq_seq_free (seq1);

  if (error != NULL)
    {
      g_printerr ("[ERROR] Writing sequence failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
    }
  g_string_free (buffer, TRUE);
  fastq_iter_free (iter1);
  if (!data.one_input)
    fastq_iter_free (iter2);

  cleanup (&data);

  return 0;
}

static void
parse_args (CallbackData   *data,
            int            *argc,
            char         ***argv)
{
  GOptionEntry entries[] =
    {
      {"output1", 'o', 0, G_OPTION_ARG_FILENAME, &data->output_path1,  "Output path for first in pair (or both if -p is not supplied)",  NULL},
      {"output2", 'p', 0, G_OPTION_ARG_FILENAME, &data->output_path2,  "Output path for second in pair (if different from -o)",          NULL},
      {"single",  's', 0, G_OPTION_ARG_FILENAME, &data->output_single, "File for single reads",                                          NULL},
      {"minsize", 'm', 0, G_OPTION_ARG_INT,      &data->min_size,      "Min read size",                                                  NULL},
      {"append",  'a', 0, G_OPTION_ARG_NONE,     &data->append,        "Append to output",                                               NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->output_path1      = NULL;
  data->output_path2      = NULL;
  data->output_single     = NULL;
  data->output_channel1   = NULL;
  data->output_channel2   = NULL;
  data->single_channel    = NULL;
  data->append            = 0;
  data->use_stdout1       = 0;
  data->use_stdout2       = 0;
  data->use_stdout_single = 0;
  data->min_size          = 0;
  data->one_input         = 0;
  data->one_output        = 0;

  context = g_option_context_new ("FILE1 [FILE2] - Removes pairs of empty sequences and separates singletons");
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
  data->input_path1 = (*argv)[1];
  if (*argc == 2)
    data->one_input = 1;
  else
    {
      data->input_path2 = (*argv)[2];
      if (!strcmp (data->input_path1, data->input_path2))
        data->one_input = 1;
    }

  if (!data->output_path1 && !data->output_path2)
    {
      data->output_path1 = g_strdup ("-");
      data->one_output   = 1;
    }
  else if (!data->output_path1)
    {
      g_printerr ("[ERROR] --output2 (-p) was provided without --output1 (-o)\n");
      exit (1);
    }
  else if (!data->output_path2 || !strcmp (data->output_path1, data->output_path2))
    data->one_output = 1;

  open_output (&data->output_channel1,
               &data->use_stdout1,
               data->output_path1,
               data->append ? "a" : "w",
               data->one_output ?
               "Opening sequence output file failed" :
               "Opening sequence output file 1 failed");

  if (!data->one_output)
    open_output (&data->output_channel2,
                 &data->use_stdout2,
                 data->output_path2,
                 data->append ? "a" : "w",
                 "Opening sequence output file 2 failed");
  else
    data->output_channel2 = data->output_channel1;

  if (data->output_single)
    open_output (&data->single_channel,
                 &data->use_stdout_single,
                 data->output_single,
                 data->append ? "a" : "w",
                 "Opening single reads output file failed");
}

static void
cleanup (CallbackData *data)
{
  close_output (&data->output_channel1,
                data->use_stdout1,
                data->one_output ?
                "Closing output file failed" :
                "Closing output file 1 failed");

  if (!data->one_output)
    close_output (&data->output_channel2,
                  data->use_stdout2,
                  "Closing output file 2 failed");

  close_output (&data->single_channel,
                data->use_stdout_single,
                "Closing single reads output file failed");

  if (data->output_path1)
    g_free (data->output_path1);
  if (data->output_path2)
    g_free (data->output_path2);
  if (data->output_single)
    g_free (data->output_single);
}

static void
open_output (GIOChannel **channel,
             int         *use_stdout,
             const char  *path,
             const char  *mode,
             const char  *error_message)
{
  GError *error = NULL;

  if (path[0] == '-' && path[1] == '\0')
    {
      *use_stdout = 1;
      *channel    = g_io_channel_unix_new (STDOUT_FILENO);
    }
  else
    {
      *channel = g_io_channel_new_file (path, mode, &error);
      if (error)
        {
          g_printerr ("[ERROR] %s: %s\n", error_message, error->message);
          exit (1);
        }
    }
  g_io_channel_set_encoding (*channel, NULL, NULL);
  g_io_channel_set_buffered (*channel, TRUE);
  g_io_channel_set_buffer_size (*channel, 4096 * 128);
}

static void
close_output (GIOChannel **channel,
              int          use_stdout,
              const char  *error_message)
{
  GError *error = NULL;

  if (*channel == NULL)
    return;
  if (!use_stdout)
    {
      g_io_channel_shutdown (*channel, TRUE, &error);
      if (error != NULL)
        {
          g_printerr ("[ERROR] %s: %s\n", error_message, error->message);
          g_error_free (error);
        }
    }
  g_io_channel_unref (*channel);
  *channel = NULL;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
