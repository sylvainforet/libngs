/**
 *
 */

#include <stdlib.h>
#include <unistd.h>

#include "ngs_fastq.h"

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  int         start;
  int         end;
  int         tot_trim;
  char       *input_path;
  char       *output_path;
  GIOChannel *out_channel;

  int         use_stdout: 1;
};

static void parse_args (CallbackData      *data,
                        int               *argc,
                        char            ***argv);

static int  iter_func  (FastqSeq          *fastq,
                        CallbackData      *data);

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
      g_printerr ("[ERROR] Iterating sequences failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
      return 1;
    }
  if (data.out_channel)
    {
      if (!data.use_stdout)
        {
          g_io_channel_shutdown (data.out_channel, TRUE, &error);
          if (error)
            {
              g_printerr ("[ERROR] Closing output file failed: %s\n", error->message);
              g_error_free (error);
              error = NULL;
            }
        }
      g_io_channel_unref (data.out_channel);
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
      {"start", 's', 0, G_OPTION_ARG_INT, &data->start, "Number of positions to trim at the start", NULL},
      {"end",   'e', 0, G_OPTION_ARG_INT, &data->end,   "Number of positions to trim at the end", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->start       = 0;
  data->end         = 0;
  data->output_path = NULL;
  data->out_channel = NULL;
  data->use_stdout  = 1;

  context = g_option_context_new ("FILE - trims the beginning and end of fastq reads");
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

  if (!data->output_path || (data->output_path[0] == '-' && data->output_path[1] == '\0'))
    data->out_channel = g_io_channel_unix_new (STDOUT_FILENO);
  else
    {
      data->out_channel = g_io_channel_new_file (data->output_path, "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening output file failed: %s\n", error->message);
          exit (1);
        }
      data->use_stdout = 0;
    }
  data->tot_trim = data->start + data->end;
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  GError  *error = NULL;
  GString *buffer;
  int      ret = 1;

  if (data->tot_trim > fastq->size)
    {
      g_printerr ("[ERROR] trimming more than sequence length (trim: %d - length: %d)\n",
                  data->tot_trim,
                  fastq->size);
      exit (1);
    }

  buffer = g_string_sized_new (512);
  buffer = g_string_append_c (buffer, '@');
  buffer = g_string_append (buffer, fastq->name + 1);
  buffer = g_string_append_c (buffer, '\n');

  fastq->seq[fastq->size - data->end] = '\0';
  buffer = g_string_append (buffer, fastq->seq + data->start);
  buffer = g_string_append_c (buffer, '\n');

  buffer = g_string_append_c (buffer, '+');
  buffer = g_string_append (buffer, fastq->name + 1);
  buffer = g_string_append_c (buffer, '\n');

  fastq->qual[fastq->size - data->end] = '\0';
  buffer = g_string_append (buffer, fastq->qual + data->start);
  buffer = g_string_append_c (buffer, '\n');

  g_io_channel_write_chars (data->out_channel,
                            buffer->str,
                            -1,
                            NULL,
                            &error);
  if (error)
    {
      g_printerr ("[ERROR] Writing trimmed sequence failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
      ret   = 0;
    }
  g_string_free (buffer, TRUE);

  return ret;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
