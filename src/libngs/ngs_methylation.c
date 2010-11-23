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
#include <unistd.h>

#include "ngs_methylation.h"


RefMethCounts*
ref_meth_counts_create (SeqDB *ref)
{
  RefMethCounts *counts;
  unsigned long  n_cg;
  unsigned long  i;

  counts = g_slice_new0 (RefMethCounts);

  n_cg = 0;
  for (i = 0; i < ref->total_size; i++)
    if (ref->seqs[i] == 'C' || ref->seqs[i] == 'G')
      ++n_cg;

  counts->meth_index = g_malloc0 (ref->total_size * sizeof (*counts->meth_index) );
  counts->meth_data  = g_malloc0 (n_cg * sizeof (*counts->meth_data) );

  n_cg = 0;
  for (i = 0; i < ref->total_size; i++)
    if (ref->seqs[i] == 'C' || ref->seqs[i] == 'G')
      counts->meth_index[i] = &counts->meth_data[n_cg++];

  return counts;
}

void
ref_meth_counts_destroy (RefMethCounts *counts)
{
  if (counts)
    {
      if (counts->meth_data)
        g_free (counts->meth_data);
      if (counts->meth_index)
        g_free (counts->meth_index);
      g_slice_free (RefMethCounts, counts);
    }
}


void
ref_meth_counts_add_path (RefMethCounts *counts,
                          SeqDB         *ref,
                          const char    *path,
                          GError       **error)
{
  GIOChannel    *channel;
  SeqDBElement  *elem      = NULL;
  GError        *tmp_error = NULL;
  char          *buffer    = NULL;
  gsize          length    = 0;
  gsize          j         = 0;
  gsize          i;

  /* Open */
  channel = g_io_channel_new_file (path, "r", &tmp_error);
  if (tmp_error)
    {
      g_propagate_error (error, tmp_error);
      return;
    }

  g_io_channel_read_to_end (channel, &buffer, &length, &tmp_error);
  if (tmp_error)
    {
      g_propagate_error (error, tmp_error);
      if (buffer)
        g_free (buffer);
      return;
    }

  for (i = 0; i < length; )
    {
      if (buffer[i] == '>')
        {
          for (j = i + 1; j < length && buffer[j] != '\n'; j++);
          buffer[j] = '\0';
          elem = g_hash_table_lookup (ref->index, buffer + i + 1);
          if (!elem)
            g_printerr ("[WARNING] Reference `%s' not found\n", buffer + i + 1);
        }
      else if (elem)
        {
          unsigned long offset;
          gsize         starts[3] = {0, 0, 0};
          int           field_idx = 0;

          for (j = i; j < length && buffer[j] != '\n'; j++)
            {
              if (buffer[j] == '\t')
                {
                  buffer[j] = '\0';
                  if (field_idx < 3)
                    starts[field_idx] = j + 1;
                  ++field_idx;
                }
            }
          buffer[j] = '\0';
          if (field_idx == 2)
            {
              offset = g_ascii_strtoll (buffer + i, NULL, 10);
              counts->meth_index[elem->offset + offset]->n_meth   += g_ascii_strtoll (buffer + starts[0], NULL, 10);
              counts->meth_index[elem->offset + offset]->n_unmeth += g_ascii_strtoll (buffer + starts[1], NULL, 10);
            }
          else if (field_idx == 3)
            {
              offset = g_ascii_strtoll (buffer + starts[0], NULL, 10);
              counts->meth_index[elem->offset + offset]->n_meth   += g_ascii_strtoll (buffer + starts[1], NULL, 10);
              counts->meth_index[elem->offset + offset]->n_unmeth += g_ascii_strtoll (buffer + starts[2], NULL, 10);
            }
          else if (field_idx != 0)
            g_printerr ("[WARNING] Could not parse meth count line\n");
        }
      i = j + 1;
    }
  if (buffer)
    g_free (buffer);

  /* Close */
  g_io_channel_shutdown (channel, TRUE, &tmp_error);
  if (tmp_error)
    {
      g_printerr ("[WARNING] Closing meth count file `%s' failed: %s\n",
                  path,
                  tmp_error->message);
      g_error_free (tmp_error);
    }
  g_io_channel_unref (channel);
}

RefMethCounts*
ref_meth_counts_load (SeqDB        *ref,
                      const char   *path,
                      GError      **error)
{
  RefMethCounts *counts;
  GError        *tmp_error = NULL;

  counts = ref_meth_counts_create (ref);
  ref_meth_counts_add_path (counts, ref, path, &tmp_error);
  if (tmp_error)
    {
      ref_meth_counts_destroy (counts);
      counts = NULL;
      g_propagate_error (error, tmp_error);
    }

  return counts;
}


void
ref_meth_counts_write (RefMethCounts *counts,
                       SeqDB         *ref,
                       const char    *path,
                       int            print_letter,
                       int            print_all,
                       GError       **error)
{
  GHashTableIter iter;
  GIOChannel    *channel;
  SeqDBElement  *elem;
  GError        *tmp_error  = NULL;
  int            use_stdout = 1;

  /* Open */
  if (!path || !*path || (path[0] == '-' && path[1] == '\0'))
    channel = g_io_channel_unix_new (STDOUT_FILENO);
  else
    {
      use_stdout = 0;
      channel = g_io_channel_new_file (path, "w", &tmp_error);
      if (tmp_error)
        {
          g_propagate_error (error, tmp_error);
          return;
        }
    }

  /* TODO could use ref_meth_counts_write_segment here */
  g_hash_table_iter_init (&iter, ref->index);
  while (g_hash_table_iter_next (&iter, NULL, (gpointer*)&elem))
    {
      GString      *buffer;
      MethCount   **ref_meth;
      unsigned long i;

      ref_meth = counts->meth_index + elem->offset;
      buffer   = g_string_new (NULL);
      g_string_printf (buffer, ">%s\n", elem->name);
      for (i = 0; i < elem->size; i++)
        {
        if (ref_meth[i])
          {
            if (print_letter)
              g_string_append_printf (buffer,
                                      "%c\t%ld\t%u\t%u\n",
                                      ref->seqs[elem->offset + i],
                                      i,
                                      ref_meth[i]->n_meth,
                                      ref_meth[i]->n_unmeth);
            else
              g_string_append_printf (buffer,
                                      "%ld\t%u\t%u\n",
                                      i,
                                      ref_meth[i]->n_meth,
                                      ref_meth[i]->n_unmeth);
          }
        else if (print_all)
          g_string_append_printf (buffer,
                                  "%c\n",
                                  ref->seqs[elem->offset + i]);
        }
      g_io_channel_write_chars (channel,
                                buffer->str,
                                buffer->len,
                                NULL,
                                &tmp_error);
      g_string_free (buffer, TRUE);

      if (tmp_error)
        {
          g_propagate_error (error, tmp_error);
          tmp_error = NULL;
          break;
        }
    }

  /* Close */
  if (!use_stdout)
    {
      g_io_channel_shutdown (channel, TRUE, &tmp_error);
      if (tmp_error)
        {
          g_printerr ("[ERROR] Closing output file `%s' failed: %s\n",
                      path,
                      tmp_error->message);
          g_error_free (tmp_error);
        }
    }
  g_io_channel_unref (channel);
}

void
ref_meth_counts_write_segment (RefMethCounts *counts,
                               SeqDB         *ref,
                               GIOChannel    *channel,
                               const char    *name,
                               unsigned long  from,
                               unsigned long  to,
                               int            print_letter,
                               int            print_all,
                               GError       **error)
{
  GError         *tmp_error = NULL;
  SeqDBElement   *elem;
  GString        *buffer;
  MethCount     **ref_meth;
  unsigned long   i;

  elem = (SeqDBElement*)g_hash_table_lookup (ref->index, name);
  if (!elem)
    return;
  if (from >= elem->size)
    return;
  to       = MIN (to, elem->size);
  ref_meth = counts->meth_index + elem->offset;
  buffer   = g_string_new (NULL);
  g_string_printf (buffer, ">%s:%lu-%lu\n", elem->name, from, to);
  for (i = from; i < to; i++)
    {
      if (ref_meth[i])
        {
          if (print_letter)
            g_string_append_printf (buffer,
                                    "%c\t%ld\t%u\t%u\n",
                                    ref->seqs[elem->offset + i],
                                    i,
                                    ref_meth[i]->n_meth,
                                    ref_meth[i]->n_unmeth);
          else
            g_string_append_printf (buffer,
                                    "%ld\t%u\t%u\n",
                                    i,
                                    ref_meth[i]->n_meth,
                                    ref_meth[i]->n_unmeth);
        }
      else if (print_all)
        g_string_append_printf (buffer,
                                "%c\n",
                                ref->seqs[elem->offset + i]);
    }
  g_io_channel_write_chars (channel,
                            buffer->str,
                            buffer->len,
                            NULL,
                            &tmp_error);
  g_string_free (buffer, TRUE);

  if (tmp_error)
    g_propagate_error (error, tmp_error);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
