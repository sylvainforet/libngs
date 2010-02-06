/**
 *
 */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ngs_fastq.h"

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char       *input_path;
  char       *output_path;
  GIOChannel *output_channel;

  int         use_stdout: 1;
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
              g_printerr ("[ERROR] Closing output file failed: %s\n", error->message);
              g_error_free (error);
              error = NULL;
            }
        }
      g_io_channel_unref (data.output_channel);
    }
  if (data.output_path)
    g_free (data.output_path);

  return 0;
}

static void
parse_args (CallbackData   *data,
            int            *argc,
            char         ***argv)
{
  GOptionEntry entries[] =
    {
      {"output" , 'o', 0, G_OPTION_ARG_FILENAME, &data->output_path, "Output path", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->output_path = NULL;
  data->use_stdout  = 1;

  context = g_option_context_new ("FILE - Splits fastq files where both mate pairs are in the same sequence");
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
  GError  *error = NULL;
  GString *buffer;
  int      ret   = 1;

  buffer = g_string_sized_new (512);

  buffer = g_string_append_c (buffer, '@');
  buffer = g_string_append (buffer, fastq->name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append_len (buffer, fastq->seq, fastq->size / 2);
  buffer = g_string_append_c (buffer, '\n');

  buffer = g_string_append_c (buffer, '+');
  buffer = g_string_append (buffer, fastq->name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append_len (buffer, fastq->qual, fastq->size / 2);
  buffer = g_string_append_c (buffer, '\n');

  fastq->name[strlen (fastq->name) - 1] = '2';
  buffer = g_string_append_c (buffer, '@');
  buffer = g_string_append (buffer, fastq->name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append (buffer, fastq->seq + fastq->size / 2);
  buffer = g_string_append_c (buffer, '\n');

  buffer = g_string_append_c (buffer, '@');
  buffer = g_string_append (buffer, fastq->name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append (buffer, fastq->qual + fastq->size / 2);
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
