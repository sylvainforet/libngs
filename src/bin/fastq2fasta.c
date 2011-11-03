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
  char       *out_seq_path;
  char       *out_qual_path;
  GIOChannel *out_seq_channel;
  GIOChannel *out_qual_channel;

  int         do_seq;
  int         do_qual;

  int         seq_use_stdout: 1;
  int         qual_use_stdout: 1;
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

  if (data.out_seq_channel)
    {
      if (!data.seq_use_stdout)
        {
          g_io_channel_shutdown (data.out_seq_channel, TRUE, &error);
          if (error)
            {
              g_printerr ("[ERROR] Closing sequence output file failed: %s\n", error->message);
              g_error_free (error);
              error = NULL;
            }
        }
      g_io_channel_unref (data.out_seq_channel);
    }
  if (data.out_qual_channel)
    {
      if (!data.qual_use_stdout)
        {
          g_io_channel_shutdown (data.out_qual_channel, TRUE, &error);
          if (error)
            {
              g_printerr ("[ERROR] Closing quality output file failed: %s\n", error->message);
              g_error_free (error);
              error = NULL;
            }
        }
      g_io_channel_unref (data.out_qual_channel);
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
      {"sequence", 's', 0, G_OPTION_ARG_NONE,     &data->do_seq,        "Extract sequences", NULL},
      {"quality" , 'q', 0, G_OPTION_ARG_NONE,     &data->do_qual,       "Extract qualities", NULL},
      {"seqout"  , 'o', 0, G_OPTION_ARG_FILENAME, &data->out_seq_path,  "Sequences file"   , NULL},
      {"qualout" , 'u', 0, G_OPTION_ARG_FILENAME, &data->out_qual_path, "Qualities file"   , NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->do_seq           = 1;
  data->do_qual          = 0;
  data->out_seq_path     = "-";
  data->out_qual_path    = NULL;
  data->out_seq_channel  = NULL;
  data->out_qual_channel = NULL;
  data->seq_use_stdout   = 1;
  data->qual_use_stdout  = 0;

  context = g_option_context_new ("FILE - Converts a fastq file to a fasta file");
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

  /* Give the filename implies do it */
  if (!data->do_seq && data->out_seq_path)
    data->do_seq = 1;
  if (!data->do_qual && data->out_qual_path)
        data->do_qual = 1;

  /* Set default output to "-" (stdout) */
  if (data->do_seq && !data->out_seq_path)
        data->out_seq_path = "-";
  if (data->do_qual && !data->out_qual_path)
        data->out_qual_path = "-";

  if (data->do_seq)
    {
      if (data->out_seq_path[0] == '-' && data->out_seq_path[1] == '\0')
        {
          data->seq_use_stdout  = 1;
          data->out_seq_channel = g_io_channel_unix_new (STDOUT_FILENO);
        }
      else
        {
          data->seq_use_stdout  = 0;
          data->out_seq_channel = g_io_channel_new_file (data->out_seq_path, "w", &error);
          if (error)
            {
              g_printerr ("[ERROR] Opening sequence output file failed: %s\n", error->message);
              exit (1);
            }
        }
    }
  if (data->do_qual)
    {
      if (data->out_qual_path[0] == '-' && data->out_qual_path[1] == '\0')
        {
          data->qual_use_stdout  = 1;
          data->out_qual_channel = g_io_channel_unix_new (STDOUT_FILENO);
        }
      else
        {
          data->qual_use_stdout  = 0;
          data->out_qual_channel = g_io_channel_new_file (data->out_qual_path, "w", &error);
          if (error)
            {
              g_printerr ("[ERROR] Opening quality output file failed: %s\n", error->message);
              exit (1);
            }
        }
    }
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  GError *error = NULL;
  int     ret   = 1;

  if (data->do_seq)
    {
      GString *buffer;

      buffer = g_string_sized_new (512);
      buffer = g_string_append_c (buffer, '>');
      buffer = g_string_append (buffer, fastq->name);
      buffer = g_string_append_c (buffer, '\n');
      buffer = g_string_append_len (buffer, fastq->seq, fastq->size);
      buffer = g_string_append_c (buffer, '\n');

      g_io_channel_write_chars (data->out_seq_channel,
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
    }
  if (data->do_qual)
    {
      GString *buffer;
      int      i;

      buffer = g_string_sized_new (512);
      buffer = g_string_append_c (buffer, '>');
      buffer = g_string_append (buffer, fastq->name);
      buffer = g_string_append_c (buffer, '\n');
      for (i = 0; i < fastq->size; i++)
        buffer = g_string_append (buffer, fastq_qual_char_2_string[(int)fastq->qual[i]]);
      buffer = g_string_append_c (buffer, '\n');

      g_io_channel_write_chars (data->out_qual_channel,
                                buffer->str,
                                -1,
                                NULL,
                                &error);
      if (error)
        {
          g_printerr ("[ERROR] Writing quality failed: %s\n", error->message);
          g_error_free (error);
          error = NULL;
          ret   = 0;
        }
      g_string_free (buffer, TRUE);
    }

  return ret;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
