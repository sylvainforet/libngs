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

int main(int argc, char **argv)
{
  return 0;
}

#if 0

#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "ngs_binseq.h"

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char              *input_path;
  GIOChannel        *input_channel;
  char              *cmd;

  unsigned int       k;
  unsigned int       k_bytes;

  int                verbose;
};

static void parse_args    (CallbackData      *data,
                           int               *argc,
                           char            ***argv);

static void kmer_bin_load (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;

  parse_args (&data, &argc, &argv);
  kmer_bin_load (&data);

  return 0;
}

static void
parse_args (CallbackData      *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      {"kmer",     'k', 0, G_OPTION_ARG_INT,      &data->k,           "K-mer size", NULL},
      {"verbose",  'v', 0, G_OPTION_ARG_NONE,     &data->verbose,     "Verbose output", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->k        = 0;
  data->verbose  = 0;

  context = g_option_context_new ("CMD FILE - Utility to manipulate kmers count files");
  g_option_context_add_main_entries (context, entries, NULL);
  if (!g_option_context_parse (context, argc, argv, &error))
    {
      g_printerr ("[ERROR] Option parsing failed: %s\n", error->message);
      exit (1);
    }
  g_option_context_free (context);

  if (*argc < 2)
    {
      g_printerr ("[ERROR] No command provided\n");
      exit (1);
    }
  data->cmd = (*argv)[1];
  if (*argc < 3)
    {
      g_printerr ("[ERROR] No input file provided\n");
      exit (1);
    }
  data->input_path = (*argv)[2];
  if (strcmp (data->cmd, "print") == 0)
    {
      if (data->verbose)
        g_printerr ("Calling the `print' command\n");
    }
  else
    {
      g_printerr ("[ERROR] Unknown command: %s\n", data->cmd);
      exit (1);
    }

  if (data->k < 1)
    {
      g_printerr ("[ERROR] The kmer size must be larger than 0\n");
      exit (1);
    }
  data->k_bytes = (data->k + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;
}

static void
kmer_bin_load (CallbackData *data)
{
  GError        *error     = NULL;
  int            use_stdin = 1;
  char          *tmp_str;
  unsigned char *buffer;
  gsize          buffer_size;

  /* Open */
  if (!data->input_path || !*data->input_path || (data->input_path[0] == '-' && data->input_path[1] == '\0'))
    data->input_channel = g_io_channel_unix_new (STDOUT_FILENO);
  else
    {
      use_stdin           = 0;
      data->input_channel = g_io_channel_new_file (data->input_path, "r", &error);
      if (error)
        {
          g_printerr ("[ERROR] failed to open input file `%s': %s\n",
                      data->input_path,
                      error->message);
          exit (1);
        }
    }
  g_io_channel_set_encoding (data->input_channel, NULL, NULL);
  /* TODO set channel buffering? */

  buffer_size      = data->k_bytes + sizeof (unsigned long int);
  buffer           = g_alloca (buffer_size);
  tmp_str          = g_alloca (data->k + 1);
  tmp_str[data->k] = '\0';

  while (1)
    {
      unsigned long int val_be;
      unsigned long int val;
      gsize             nb_read;
      GIOStatus         status;

      status = g_io_channel_read_chars (data->input_channel,
                                        (char*)buffer,
                                        buffer_size,
                                        &nb_read,
                                        &error);
      if (error)
        {
          g_printerr ("[ERROR] Reading from input file `%s' failed: %s\n",
                      data->input_path,
                      error->message);
          exit (1);
        }
      if (status == G_IO_STATUS_EOF || status != G_IO_STATUS_NORMAL)
        break;
      if (nb_read != buffer_size)
        {
          g_printerr ("[ERROR] Expected %ld bytes, read only %ld\n",
                      buffer_size,
                      nb_read);
          exit (1);
        }
      bin_to_char_prealloc (tmp_str, buffer, data->k);
      val_be = *(unsigned long int*)(buffer + data->k_bytes);
      val    = GULONG_FROM_BE (val_be);
      g_print ("%s %ld\n", tmp_str, val);

    }

  /* Close */
  if (!use_stdin)
    {
      g_io_channel_shutdown (data->input_channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing input file `%s' failed: %s\n",
                      data->input_path,
                      error->message);
          g_error_free (error);
        }
    }
  g_io_channel_unref (data->input_channel);
}

#endif

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
