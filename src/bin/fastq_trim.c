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
  int         keep;

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
      {"out",   'o', 0, G_OPTION_ARG_FILENAME, &data->output_path, "Output file", NULL},
      {"start", 's', 0, G_OPTION_ARG_INT,  &data->start, "Number of positions to trim at the start", NULL},
      {"end",   'e', 0, G_OPTION_ARG_INT,  &data->end,   "Number of positions to trim at the end", NULL},
      {"qual",  'q', 0, G_OPTION_ARG_INT,  &data->qual,  "Minimum quality", NULL},
      {"len",   'l', 0, G_OPTION_ARG_INT,  &data->len,   "Minimum length", NULL},
      {"keep",  'k', 0, G_OPTION_ARG_NONE, &data->keep,  "Keep an pseudo-entry for too small reads", NULL},
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
  data->keep           = 0;

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
  data->tot_trim = data->start + data->end;
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  int end;

  if (data->tot_trim > fastq->size)
    {
      g_printerr ("[ERROR] trimming more than sequence length "
                  "(trim: %d - length: %d)\n",
                  data->tot_trim, fastq->size);
      exit (1);
    }

  end = fastq->size - data->end;
  if (data->qual > 0)
    {
      int i;
      for (i = 0; i < fastq->size; i++)
        if (fastq->qual[i] - FASTQ_QUAL_0 < data->qual)
          break;
      if (i < end)
        end = i;
    }

  if (end <= data->start ||
      end - data->start < data->len)
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
                      data->start,
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
