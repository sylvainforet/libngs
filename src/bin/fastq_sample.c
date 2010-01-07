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
  GRand      *rand;
  FastqSeq  **samples;

  int         n;
  int         t;

  int         use_stdout: 1;
};

static int  iter_func     (FastqSeq       *fastq,
                           CallbackData   *data);

static void parse_args    (CallbackData   *data,
                           int            *argc,
                           char         ***argv);

static void print_samples (CallbackData   *data);

static void free_samples  (CallbackData   *data);

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

  print_samples (&data);

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
  g_rand_free (data.rand);
  free_samples (&data);

  return 0;
}

static void
parse_args (CallbackData   *data,
            int            *argc,
            char         ***argv)
{
  GOptionEntry entries[] =
    {
      {"number", 'n', 0, G_OPTION_ARG_INT, &data->n, "Number of samples", NULL},
      {"out",    'o', 0, G_OPTION_ARG_FILENAME, &data->output_path, "Output file", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->n              = 1;
  data->t              = 0;
  data->output_path    = "-";
  data->output_channel = NULL;
  data->samples        = NULL;
  data->use_stdout     = 1;

  context = g_option_context_new ("FILE - Samples reads from a fastq file");
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

  if (data->output_path[0] == '-' && data->output_path[1] == '\0')
    data->output_channel = g_io_channel_unix_new (STDOUT_FILENO);
  else
    {
      data->use_stdout     = 0;
      data->output_channel = g_io_channel_new_file (data->output_path, "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening quality output file failed: %s\n", error->message);
          exit (1);
        }
    }

  data->samples = g_malloc (data->n * sizeof (*data->samples));
  data->rand    = g_rand_new ();
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  if (data->t < data->n)
    data->samples[data->t] = fastq_seq_copy (fastq);
  else
    {
      const int rnd = g_rand_int_range (data->rand,
                                        0, data->t);

      if (rnd < data->n)
        data->samples[rnd] = fastq_seq_copy (fastq);
    }
  data->t++;

  return 1;
}

static void
print_samples (CallbackData *data)
{
  int i;

  if (data->t < data->n)
    {
      g_printerr ("[ERROR] Found less records than requested samples "
                  "(found %d, requested %d)\n",
                  data->t, data->n);
      return;
    }
  for (i = 0; i < data->n; i++)
    {
      GError  *error = NULL;
      GString *buffer;

      buffer = g_string_sized_new (512);
      buffer = g_string_append_c (buffer, '@');
      buffer = g_string_append (buffer, data->samples[i]->name);
      buffer = g_string_append_c (buffer, '\n');

      buffer = g_string_append (buffer, data->samples[i]->seq);
      buffer = g_string_append_c (buffer, '\n');

      buffer = g_string_append_c (buffer, '+');
      buffer = g_string_append (buffer, data->samples[i]->name);
      buffer = g_string_append_c (buffer, '\n');

      buffer = g_string_append (buffer, data->samples[i]->qual);
      buffer = g_string_append_c (buffer, '\n');

      g_io_channel_write_chars (data->output_channel,
                                buffer->str,
                                -1,
                                NULL,
                                &error);
      g_string_free (buffer, TRUE);
      if (error)
        {
          g_printerr ("[ERROR] Writing sampled sequence failed: %s\n",
                      error->message);
          g_error_free (error);
          return;
        }
    }
}

static void
free_samples (CallbackData *data)
{
  int n;
  int i;

  if (data->samples)
    {
      n = MIN (data->t, data->n);
      for (i = 0; i < n; i++)
        fastq_seq_free (data->samples[i]);
      g_free (data->samples);
    }
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
