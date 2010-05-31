/**
 *
 */

#include <stdlib.h>
#include <unistd.h>

#include "ngs_binseq.h"

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char              *input_path;
  GIOChannel        *input_channel;

  unsigned char     *tmp_word;
  char              *tmp_str;

  unsigned int       k;
  unsigned int       k_bytes;
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

  if (data.tmp_word)
    g_free (data.tmp_word);
  if (data.tmp_str)
    g_free (data.tmp_str);

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
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->k        = 0;
  data->tmp_word = NULL;
  data->tmp_str  = NULL;

  context = g_option_context_new ("FILE - blah blah fucking blah");
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

  if (data->k < 1)
    {
      g_printerr ("[ERROR] The kmer size must be larger than 0\n");
      exit (1);
    }

  data->k_bytes          = (data->k + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;
  data->tmp_word         = g_malloc (data->k_bytes);
  data->tmp_str          = g_malloc (data->k + 1);
  data->tmp_str[data->k] = '\0';

  data->input_path = (*argv)[1];
}

static void
kmer_bin_load (CallbackData *data)
{
  GError           *error     = NULL;
  int               use_stdin = 1;
  unsigned long int val_be;
  unsigned long int val;

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

  while (1)
    {
      GIOStatus status;

      status = g_io_channel_read_chars (data->input_channel,
                                        (char*)data->tmp_word,
                                        data->k_bytes,
                                        NULL,
                                        &error);
      if (error)
        {
          g_printerr ("[ERROR] Reading from input file `%s' failed: %s\n",
                      data->input_path,
                      error->message);
          g_error_free (error);
          exit (1);
        }
      if (status != G_IO_STATUS_NORMAL)
        break;
      status = g_io_channel_read_chars (data->input_channel,
                                        (char*)&val_be,
                                        sizeof (unsigned long int),
                                        NULL,
                                        &error);
      if (error)
        {
          g_printerr ("[ERROR] Reading from input file `%s' failed: %s\n",
                      data->input_path,
                      error->message);
          g_error_free (error);
          exit (1);
        }
      val = GULONG_FROM_BE (val_be);
      bin_to_char_prealloc (data->tmp_str, data->tmp_word, data->k);
      g_print ("%s %ld\n", data->tmp_str, val);
      if (status != G_IO_STATUS_NORMAL)
        break;
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

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
