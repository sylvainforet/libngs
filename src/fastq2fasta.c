/**
 *
 */

#include <stdlib.h>
#include <unistd.h>

#include "fastq.h"

typedef struct _Fastq2FastaData Fastq2FastaData;

struct _Fastq2FastaData
{
  char       *input_path;
  char       *out_seq_path;
  char       *out_qual_path;
  GIOChannel *out_seq_channel;
  GIOChannel *out_qual_channel;
  int         do_seq;
  int         do_qual;
};

static int  iter_func  (FastqSeq          *fastq,
                        Fastq2FastaData   *data);

static void parse_args (Fastq2FastaData   *data,
                        int               *argc,
                        char            ***argv);

int
main (int    argc,
      char **argv)
{
  Fastq2FastaData data;
  GError         *error = NULL;

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
      g_io_channel_shutdown (data.out_seq_channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing sequence output file failed: %s\n", error->message);
          g_error_free (error);
          error = NULL;
        }
      g_io_channel_unref (data.out_seq_channel);
    }
  if (data.out_qual_channel)
    {
      g_io_channel_shutdown (data.out_qual_channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing quality output file failed: %s\n", error->message);
          g_error_free (error);
          error = NULL;
        }
      g_io_channel_unref (data.out_qual_channel);
    }

  return 0;
}

static void
parse_args (Fastq2FastaData   *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      {"sequence", 's', 0, G_OPTION_ARG_NONE,   &data->do_seq,        "Extract sequences", NULL},
      {"quality" , 'q', 0, G_OPTION_ARG_NONE,   &data->do_qual,       "Extract qualities", NULL},
      {"seqout"  , 'o', 0, G_OPTION_ARG_STRING, &data->out_seq_path,  "Sequences file"   , NULL},
      {"qualout" , 'u', 0, G_OPTION_ARG_STRING, &data->out_qual_path, "Qualities file"   , NULL},
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

  context = g_option_context_new ("FILE - Converts a fastq file to a fasta file");
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

  if (!data->do_seq && data->out_seq_path)
    data->do_seq = 1;
  if (!data->do_qual && data->out_qual_path)
        data->do_qual = 1;
  if (data->do_seq && !data->out_seq_path)
        data->out_seq_path = "-";
  if (data->do_qual && !data->out_qual_path)
        data->out_qual_path = "-";
  if (data->do_seq)
    {
      if (data->out_seq_path[0] == '-' && data->out_seq_path[1] == '\0')
        data->out_seq_channel = g_io_channel_unix_new (STDOUT_FILENO);
      else
        {
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
        data->out_qual_channel = g_io_channel_unix_new (STDOUT_FILENO);
      else
        {
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
iter_func (FastqSeq        *fastq,
           Fastq2FastaData *data)
{
  GError *error = NULL;
  int     ret   = 1;

  if (data->do_seq)
    {
      char *buffer;

      buffer = g_strdup_printf (">%s\n%s\n",
                                fastq->name + 1,
                                fastq->seq);
      g_io_channel_write_chars (data->out_seq_channel,
                                buffer,
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
      g_free (buffer);
    }
  if (data->do_qual)
    {
      GString *buffer;
      char    *tmp;

      buffer = g_string_sized_new (512);
      buffer = g_string_append_c (buffer, '>');
      buffer = g_string_append (buffer, fastq->name + 1);
      buffer = g_string_append_c (buffer, '\n');
      tmp    = fastq->qual - 1;
      while (*++tmp)
        buffer = g_string_append (buffer, fastq_qual_char_2_string[(int)*tmp]);
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
  fastq_seq_free (fastq);

  return ret;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
