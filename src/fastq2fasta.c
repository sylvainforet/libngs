/**
 *
 */

#include <unistd.h>

#include "fastq.h"

typedef struct _Fastq2FastaData Fastq2FastaData;

struct _Fastq2FastaData
{
  GIOChannel *out_channel;
};

static int iter_func (FastqSeq        *fastq,
                      Fastq2FastaData *data);

int
main (int    argc,
      char **argv)
{
  Fastq2FastaData data;
  GError         *error = NULL;

  if (argc < 2)
    {
      g_print ("Usage: fastq2fasta FILE\n");
      return 1;
    }

  data.out_channel = g_io_channel_unix_new (STDOUT_FILENO);
  if (!data.out_channel)
    return 1;

  iter_fastq (argv[1],
              (FastqIterFunc)iter_func,
              (void*)&data,
              &error);
  if (error)
    {
      g_printerr ("[ERROR] %s\n", error->message);
      g_error_free (error);
      g_io_channel_unref (data.out_channel);
      return -1;
    }
  g_io_channel_unref (data.out_channel);

  return 0;
}

static int
iter_func (FastqSeq        *fastq,
           Fastq2FastaData *data)
{
  g_print(">%s\n%s\n",
          fastq->name + 1,
          fastq->seq);
  fastq_seq_free (fastq);

  return 1;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
