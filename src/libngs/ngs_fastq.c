/**
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "ngs_fastq.h"
#include "ngs_fastq_flex.h"
#include "ngs_utils.h"


static char *fastq_parser_name = NULL;

static void iter_fastq_simple (const char   *path,
                               FastqIterFunc func,
                               void         *data,
                               GError      **error);

static void iter_fastq_ugly   (const char   *path,
                               FastqIterFunc func,
                               void         *data,
                               GError      **error);

static int  iter_load_db      (FastqSeq     *fastq,
                               FastqDB      *db);

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

FastqSeq*
fastq_seq_copy (FastqSeq *fastq)
{
  FastqSeq *seq = NULL;

  if (fastq)
    {
      seq       = fastq_seq_new ();
      seq->name = strdup (fastq->name);
      seq->seq  = strdup (fastq->seq);
      seq->qual = strdup (fastq->qual);
      seq->size = fastq->size;
    }

  return seq;
}

void
iter_fastq (const char   *path,
            FastqIterFunc func,
            void         *data,
            GError      **error)
{
  if (fastq_parser_name == NULL || g_strcmp0 (fastq_parser_name, "flex") == 0)
    iter_fastq_flex (path, func, data, error);
  else if (g_strcmp0 (fastq_parser_name, "simple") == 0)
    iter_fastq_simple (path, func, data, error);
  else if (g_strcmp0 (fastq_parser_name, "ugly") == 0)
    iter_fastq_ugly (path, func, data, error);
  else if (g_strcmp0 (fastq_parser_name, "flexugly") == 0)
    iter_fastq_flex_ugly (path, func, data, error);
  else
    {
      g_set_error (error,
                   NGS_ERROR,
                   NGS_UNKNOWN_ERROR,
                   "Unknown fastq parser: `%s'",
                   fastq_parser_name);
    }
}

static void
iter_fastq_simple (const char   *path,
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
  int         ret;

  if (path[0] == '-' && path[1] == '\0')
    {
      channel = g_io_channel_unix_new (STDIN_FILENO);
      if (!channel)
        {
          g_set_error (error,
                       NGS_ERROR,
                       NGS_IO_ERROR,
                       "Could not open stdin for reading");
          return;
        }
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
          fastq->name = strdup (line + 1);
          g_free (line);

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

          /* Size */
          fastq->size = MIN (strlen (fastq->seq),
                             strlen (fastq->qual));

          /* Callback */
          ret = func (fastq, data);
          fastq_seq_free (fastq);
          if (!ret)
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

static void
iter_fastq_ugly  (const char    *path,
                   FastqIterFunc func,
                   void         *data,
                   GError      **error)
{
#define BUFFER_SIZE (1024 * 16)
#define MAX_STRING_SIZE 256
  char      buffer[BUFFER_SIZE];
  char      name[MAX_STRING_SIZE + 1];
  char      seq[MAX_STRING_SIZE + 1];
  char      qual[MAX_STRING_SIZE + 1];
  FILE     *file;
  FastqSeq  current;
  FastqSeq *fastq    = NULL;
  int       name_idx = 0;
  int       seq_idx  = 0;
  int       qual_idx = 0;
  int       status   = 0;
  int       i;

  if (path[0] == '-' && path[1] == '\0')
    file = stdin;
  else
    file = fopen (path, "r");
  if (file == NULL)
    {
      g_set_error (error,
                   NGS_ERROR,
                   NGS_IO_ERROR,
                   "Could not open input file `%s'",
                   path);
      return;
    }

  do
    {
      const size_t bytes_read = fread (buffer, sizeof (char), BUFFER_SIZE, file);

      for (i = 0; i < bytes_read; i++)
        {
          if (status == 0 && buffer[i] == '@')
            {
              if (fastq)
                {
                  int ret;

                  fastq->size = MIN (seq_idx, qual_idx);
                  ret         = func (fastq, data);
                  if (!ret)
                    goto cleanup;
                }
              else
                {
                  fastq       = &current;
                  fastq->name = name;
                  fastq->seq  = seq;
                  fastq->qual = qual;
                }
              name_idx = 0;
              seq_idx  = 0;
              qual_idx = 0;
              status   = 1;
            }
          /* `@' before name */
          if (status == 1)
            status = 2;
          /* name */
          else if (status == 2)
            {
              if (buffer[i] == '\n')
                {
                  name[name_idx] = '\0';
                  status = 3;
                }
              else if (name_idx == MAX_STRING_SIZE)
                {
                  if (buffer[i] == '\n')
                    {
                      name[name_idx] = '\0';
                      status = 3;
                    }
                }
              else
                {
                  name[name_idx] = buffer[i];
                  ++name_idx;
                }
            }
          /* sequence */
          else if (status == 3)
            {
              if (buffer[i] == '\n')
                {
                  seq[seq_idx] = '\0';
                  status = 4;
                }
              else if (seq_idx == MAX_STRING_SIZE)
                {
                  if (buffer[i] == '\n')
                    {
                      seq[seq_idx] = '\0';
                      status = 4;
                    }
                }
              else
                {
                  seq[seq_idx] = buffer[i];
                  ++seq_idx;
                }
            }
          /* name again */
          else if (status == 4)
            {
              if (buffer[i] == '\n')
                status = 5;
            }
          /* quality */
          else if (status == 5)
            {
              if (buffer[i] == '\n')
                {
                  qual[qual_idx] = '\0';
                  status = 0;
                }
              else if (qual_idx == MAX_STRING_SIZE)
                {
                  if (buffer[i] == '\n')
                    {
                      qual[qual_idx] = '\0';
                      status = 0;
                    }
                }
              else
                {
                  qual[qual_idx] = buffer[i];
                  ++qual_idx;
                }
            }
        }
      if (bytes_read < BUFFER_SIZE)
        {
          if (feof (file))
            break;
          if (ferror (file))
            {
              g_set_error (error,
                           NGS_ERROR,
                           NGS_IO_ERROR,
                           "Error while reading `%s'",
                           path);
              goto cleanup;
            }
        }
    }
  while (1);
  if ((status == 0 || status == 5) && fastq)
    {
      fastq->size = MIN (seq_idx, qual_idx);
      func (fastq, data);
    }

cleanup:
  fclose (file);

#undef BUFFER_SIZE
#undef MAX_NAME_SIZE
}

GOptionGroup*
get_fastq_option_group (void)
{
  GOptionEntry entries[] =
    {
      {"fastq_parser_name", 0, 0, G_OPTION_ARG_STRING, &fastq_parser_name, "Name of the parser", NULL},
      {NULL}
    };
  GOptionGroup *option_group;

  option_group = g_option_group_new ("fastq",
                                     "fastq parser options",
                                     "Show fastq parser options",
                                     NULL,
                                     NULL);
  g_option_group_add_entries (option_group, entries);

  return option_group;
}

FastqDB*
fastq_db_new (void)
{
  FastqDB *db;

  db        = g_malloc (sizeof (*db));
  db->index = g_hash_table_new_full (g_str_hash,
                                     g_str_equal,
                                     g_free,
                                     NULL);
  db->seqs       = NULL;
  db->quals      = NULL;
  db->seq_size   = 0;
  db->n_seqs     = 0;
  db->alloc_size = 0;

  return db;
}

void
fastq_db_load (FastqDB     *db,
               const char  *path,
               GError     **error)
{
  iter_fastq (path,
              (FastqIterFunc)iter_load_db,
              db,
              error);
}

void
fastq_db_free (FastqDB *db)
{
  if (db)
    {
      if (db->index)
        g_hash_table_destroy (db->index);
      if (db->seqs)
        g_free (db->seqs);
      if (db->quals)
        g_free (db->quals);
      g_free (db);
    }
}

#define SIZE_INCREMENT (1024 * 1024)
static int
iter_load_db (FastqSeq *fastq,
              FastqDB  *db)
{
  if (db->seq_size == 0)
    db->seq_size = fastq->size;
  else if (db->seq_size != fastq->size)
    {
      g_printerr ("[WARNING] Sequence `%s' has different size, ignoring\n",
                  fastq->name);
      return 1;
    }
  db->n_seqs++;
  if (db->n_seqs * db->seq_size * sizeof (*db->seqs) > db->alloc_size)
    {
      db->alloc_size += SIZE_INCREMENT;
      db->seqs        = g_realloc (db->seqs, db->alloc_size);
      db->quals       = g_realloc (db->quals, db->alloc_size);
    }
  g_hash_table_insert (db->index,
                       strdup (fastq->name),
                       GUINT_TO_POINTER (db->n_seqs - 1));

  return 1;
}
#undef SIZE_INCREMENT

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
