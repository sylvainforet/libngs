/**
 *
 */

#include <unistd.h>

#include "fastq.h"


FastqSeq*
fastq_seq_new (void)
{
  FastqSeq *fastq;

  fastq = g_slice_new0 (FastqSeq);

  return fastq;
}

void
fastq_seq_free (FastqSeq *fastq)
{
  if (fastq)
    {
      if (fastq->name)
        g_free (fastq->name);
      if (fastq->seq)
        g_free (fastq->seq);
      if (fastq->qual)
        g_free (fastq->qual);
      g_slice_free (FastqSeq, fastq);
    }
}

void
iter_fastq (char         *path,
            FastqIterFunc func,
            void         *data,
            GError      **error)
{
  GIOChannel *channel;
  GError     *tmp_err = NULL;
  FastqSeq   *fastq   = NULL;
  char       *line    = NULL;
  gsize       length;
  gsize       endl;

  if (path[0] == '-' && path[1] == '\0')
    {
      channel = g_io_channel_unix_new (STDIN_FILENO);
      if (!channel) /* TODO raise an error here */
        return;
    }
  else
    {
      channel = g_io_channel_new_file (path, "r", &tmp_err);

      if (tmp_err != NULL)
        {
          g_propagate_error (error, tmp_err);
          if (channel)
            g_io_channel_unref (channel);
          return;
        }
    }

  while (G_IO_STATUS_NORMAL == g_io_channel_read_line (channel, &line, &length, &endl, &tmp_err))
    {
      if (*line == '@')
        {
          /* Sequence header */
          line[endl]  = '\0';
          fastq       = fastq_seq_new ();
          fastq->name = line;

          /* Sequence */
          g_io_channel_read_line (channel, &line, &length, &endl, &tmp_err);
          if (tmp_err)
            goto error;
          line[endl] = '\0';
          fastq->seq = line;

          /* Quality header */
          g_io_channel_read_line (channel, &line, &length, &endl, &tmp_err);
          if (tmp_err)
            goto error;
          g_free (line);

          /* Quality */
          g_io_channel_read_line (channel, &line, &length, &endl, &tmp_err);
          if (tmp_err)
            goto error;
          line[endl]  = '\0';
          fastq->qual = line;

          /* Callback */
          if (!func (fastq, data))
            break;
        }
    }
  goto cleanup;

error:
  if (line)
    g_free (line);
  if (fastq)
    fastq_seq_free (fastq);
  g_propagate_error (error, tmp_err);

cleanup:
  g_io_channel_unref (channel);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
