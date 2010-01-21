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
#if 0
  GIOChannel    *channel;
  SeqDBElement  *elem      = NULL;
  GError        *tmp_error = NULL;
  char          *line      = NULL;
  gsize          length;
  gsize          endl;

  /* Open */
  channel = g_io_channel_new_file (path, "r", &tmp_error);
  if (tmp_error)
    {
      g_propagate_error (error, tmp_error);
      return;
    }

  /* Parse */
  while (G_IO_STATUS_NORMAL == g_io_channel_read_line (channel, &line, &length, &endl, &tmp_error))
    {
      line[endl] = '\0';
      if (line[0] == '>')
        {
          elem = g_hash_table_lookup (ref->index, line + 1);
          if (!elem)
            g_printerr ("[WARNING] Reference `%s' not found\n", line);
        }
      else if (elem)
        {
          char        **values;
          int           n_fields;
          unsigned long offset;

          values = g_strsplit (line, "\t", -1);

          for (n_fields = 0; values[n_fields]; n_fields++);
          if (n_fields == 3)
            {
              offset = g_ascii_strtoll (values[0], NULL, 10);
              counts->meth_index[elem->offset + offset]->n_meth   += g_ascii_strtoll (values[1], NULL, 10);
              counts->meth_index[elem->offset + offset]->n_unmeth += g_ascii_strtoll (values[2], NULL, 10);
            }
          else if (n_fields == 4)
            {
              offset = g_ascii_strtoll (values[1], NULL, 10);
              counts->meth_index[elem->offset + offset]->n_meth   += g_ascii_strtoll (values[2], NULL, 10);
              counts->meth_index[elem->offset + offset]->n_unmeth += g_ascii_strtoll (values[3], NULL, 10);
            }
          else if (n_fields != 1)
            g_printerr ("[WARNING] Could not parse meth count line\n");

          g_strfreev (values);
        }
      g_free (line);
      line = NULL;
    }
  if (line)
    g_free (line);
  if (tmp_error)
    {
      g_propagate_error (error, tmp_error);
      tmp_error = NULL;
    }

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
#endif /* 0 */
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
              if (field_idx > 4)
                break;
              if (buffer[i] == '\t')
                {
                  buffer[i] = '\0';
                  ++field_idx;
                  if (field_idx < 4)
                    starts[field_idx] = j + 1;
                }
            }
          buffer[j] = '\0';
          if (field_idx == 1)
            {
              offset = g_ascii_strtoll (buffer + i, NULL, 10);
              counts->meth_index[elem->offset + offset]->n_meth   += g_ascii_strtoll (buffer + i + starts[0], NULL, 10);
              counts->meth_index[elem->offset + offset]->n_unmeth += g_ascii_strtoll (buffer + i + starts[1], NULL, 10);
            }
          else if (field_idx == 2)
            {
              offset = g_ascii_strtoll (buffer + i + starts[0], NULL, 10);
              counts->meth_index[elem->offset + offset]->n_meth   += g_ascii_strtoll (buffer + i + starts[1], NULL, 10);
              counts->meth_index[elem->offset + offset]->n_unmeth += g_ascii_strtoll (buffer + i + starts[2], NULL, 10);
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

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
