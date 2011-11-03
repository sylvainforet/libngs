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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "ngs_fastq.h"
#include "ngs_fastq_flex.h"
#include "ngs_utils.h"


char  fastq_qual0          = FASTQ_QUAL_0;
char  fastq_qual_min_str[] = FASTQ_QUAL_MIN_STR;


static char*    fastq_parser_name = NULL;

static void     iter_fastq_simple      (const char     *path,
                                        FastqIterFunc   func,
                                        void           *data,
                                        GError        **error);

static void     iter_fastq_ugly        (const char     *path,
                                        FastqIterFunc   func,
                                        void           *data,
                                        GError        **error);

static gboolean fastq_parse_qual0      (const gchar    *option_name,
                                        const gchar    *value,
                                        gpointer        data,
                                        GError        **error);

static void     init_char2string_table (char            qual0);

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
      seq->name = g_strdup (fastq->name);
      seq->seq  = g_strdup (fastq->seq);
      seq->qual = g_strdup (fastq->qual);
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
  if (fastq_parser_name == NULL || strcmp (fastq_parser_name, "flex") == 0)
    iter_fastq_flex (path, func, data, error);
  else if (strcmp (fastq_parser_name, "simple") == 0)
    iter_fastq_simple (path, func, data, error);
  else if (strcmp (fastq_parser_name, "ugly") == 0)
    iter_fastq_ugly (path, func, data, error);
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
      {"fastq_parser_name", 0, 0, G_OPTION_ARG_STRING,   &fastq_parser_name, "Name of the parser",               NULL},
      {"fastq_qual0"      , 0, 0, G_OPTION_ARG_CALLBACK, &fastq_parse_qual0, "The character encoding quality 0", NULL},
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

FastqIter*
fastq_iter_new (const char *path,
                GError    **error)
{
  FastqIter *iter;

  iter = g_slice_new (FastqIter);
  iter->private = fastq_iter_new_flex (path, error);

  return iter;
}

FastqSeq*
fastq_iter_next (FastqIter *iter)
{
  return fastq_iter_next_flex ((FastqIterFlex*)iter->private);
}

void
fastq_iter_free (FastqIter *iter)
{
  if (iter)
    {
      if (iter->private)
        fastq_iter_free_flex ((FastqIterFlex*)iter->private);
      g_slice_free (FastqIter, iter);
    }
}

void
fastq_write (GIOChannel *channel,
             GString    *buffer,
             char       *name,
             char       *seq,
             char       *qual,
             GError    **error)
{
  GError  *tmp_err     = NULL;
  int      free_buffer = 0;

  if (buffer == NULL)
    {
      buffer = g_string_sized_new (512);
      free_buffer = 1;
    }
  else
    g_string_truncate (buffer, 0);

  buffer = g_string_append_c (buffer, '@');
  buffer = g_string_append (buffer, name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append (buffer, seq);
  buffer = g_string_append (buffer, "\n+\n");
  buffer = g_string_append (buffer, qual);
  buffer = g_string_append_c (buffer, '\n');

  g_io_channel_write_chars (channel,
                            buffer->str,
                            -1,
                            NULL,
                            &tmp_err);
  if (tmp_err != NULL)
    g_propagate_error (error, tmp_err);

  if (free_buffer)
    g_string_free (buffer, TRUE);
}

void
fastq_write_fragment (GIOChannel *channel,
                      GString    *buffer,
                      char       *name,
                      char       *seq,
                      char       *qual,
                      int         start,
                      int         end,
                      GError    **error)
{
  GError  *tmp_err     = NULL;
  int      free_buffer = 0;

  if (buffer == NULL)
    {
      buffer = g_string_sized_new (512);
      free_buffer = 1;
    }
  else
    g_string_truncate (buffer, 0);

  buffer = g_string_append_c (buffer, '@');
  buffer = g_string_append (buffer, name);
  buffer = g_string_append_c (buffer, '\n');
  buffer = g_string_append_len (buffer, seq + start, end - start);
  buffer = g_string_append (buffer, "\n+\n");
  buffer = g_string_append_len (buffer, qual + start, end - start);
  buffer = g_string_append_c (buffer, '\n');

  g_io_channel_write_chars (channel,
                            buffer->str,
                            -1,
                            NULL,
                            &tmp_err);
  if (tmp_err != NULL)
    g_propagate_error (error, tmp_err);

  if (free_buffer)
    g_string_free (buffer, TRUE);
}

static gboolean
fastq_parse_qual0 (const gchar  *option_name,
                   const gchar  *value,
                   gpointer      data,
                   GError      **error)
{
  if (value == NULL || strlen (value) > 1)
    {
      g_set_error (error,
                   NGS_ERROR,
                   NGS_ARG_ERROR,
                   "Option `--fastq_qual0' expects a single ascii character (got: %s)",
                   value);
      return FALSE;
    }
  fastq_qual0 = value[0];
  switch (fastq_qual0)
    {
      case '@':
      case '%':
      case '!':
          break;
      default:
          g_printerr ("[WARNING] Unconventional zero quality character: %s\n", value);
          if (fastq_qual0 > 88)
            {
              g_set_error (error,
                           NGS_ERROR,
                           NGS_ARG_ERROR,
                           "fastq_qual0 (%s) value too high, can only encode %d quality values",
                           value,
                           128 - fastq_qual0);
              return FALSE;
            }
          break;
    }
  if (fastq_qual0 != FASTQ_QUAL_0)
    {
      init_char2string_table (fastq_qual0);
      fastq_qual_min_str[0] = fastq_qual0 + 2;
    }

  return TRUE;
}

static void
init_char2string_table (char qual0)
{
  int i;

  for (i = 0; i < 128; ++i)
    {
      const int q = i - qual0;

      fastq_qual_char_2_string[i][0] = ' ';
      if (q <= 0 || q > 99)
        {
          fastq_qual_char_2_string[i][1] = '0';
          fastq_qual_char_2_string[i][2] = ' ';
        }
      else if (q <= 9)
        {
          fastq_qual_char_2_string[i][1] = '0' + q;
          fastq_qual_char_2_string[i][2] = ' ';
        }
      else
        {
          fastq_qual_char_2_string[i][1] = '0' + q / 10;
          fastq_qual_char_2_string[i][2] = '0' + q % 10;
        }
    }
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
