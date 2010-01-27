/**
 *
 */

#include <stdlib.h>
#include <string.h>

#include "ngs_fasta.h"
#include "ngs_fasta_flex.h"
#include "ngs_utils.h"


static char *fasta_parser_name = NULL;

static void iter_load_dict (FastaSeq   *fasta,
                            GHashTable *dict);

FastaSeq*
fasta_seq_new (void)
{
  FastaSeq *fasta;

  fasta = g_slice_new0 (FastaSeq);

  return fasta;
}

void
fasta_seq_free (FastaSeq *fasta)
{
  if (fasta)
    {
      if (fasta->name)
        g_free (fasta->name);
      if (fasta->seq)
        g_free (fasta->seq);
      g_slice_free (FastaSeq, fasta);
    }
}

void
iter_fasta (const char   *path,
            FastaIterFunc func,
            void         *data,
            GError      **error)
{
  if (fasta_parser_name == NULL || strcmp (fasta_parser_name, "flex") == 0)
    iter_fasta_flex (path, func, data, error);
  else
    {
      g_set_error (error,
                   NGS_ERROR,
                   NGS_UNKNOWN_ERROR,
                   "Unknown fasta parser: `%s'",
                   fasta_parser_name);
    }
}

GHashTable*
load_fasta_dict (char *path)
{
  GError     *error = NULL;
  GHashTable *data;

  data = g_hash_table_new_full (g_str_hash,
                                g_str_equal,
                                g_free,
                                g_free);
  iter_fasta (path,
              (FastaIterFunc)iter_load_dict,
              data,
              &error);
  if (error)
    {
      g_printerr ("[ERROR] Problem while loading fasta dict: %s\n", error->message);
      g_hash_table_destroy (data);
      return NULL;
    }

  return data;
}

static void
iter_load_dict (FastaSeq   *fasta,
                GHashTable *dict)
{
  g_hash_table_insert (dict,
                       strdup (fasta->name),
                       strdup (fasta->seq));
}

GOptionGroup*
get_fasta_option_group (void)
{
  GOptionEntry entries[] =
    {
      {"fasta_parser_name", 0, 0, G_OPTION_ARG_STRING, &fasta_parser_name, "Name of the fasta parser", NULL},
      {NULL}
    };
  GOptionGroup *option_group;

  option_group = g_option_group_new ("fasta",
                                     "fasta parser options",
                                     "Show fasta parser options",
                                     NULL,
                                     NULL);
  g_option_group_add_entries (option_group, entries);

  return option_group;
}

char*
fasta_write_to_buffer (FastaSeq     *seq,
                       unsigned int line_size)
{
  GString      *buffer;
  unsigned long size;
  unsigned long i;

  size = 3                     + /* for the '>', the first '\n' and the last '\0' */
         strlen (seq->name)    + /* the name */
         seq->size             + /* the sequence size */
         (seq->size + line_size - 1) / line_size; /* the '\n's */

  buffer = g_string_sized_new (size);
  buffer = g_string_append_c (buffer, '>');
  buffer = g_string_append (buffer, seq->name);
  i      = 0;
  do
    {
      buffer = g_string_append_c (buffer, '\n');
      buffer = g_string_append_len (buffer, seq->seq + i, line_size);
      i     += line_size;
    }
  while (i < seq->size);
  buffer = g_string_append_c (buffer, '\n');

  return g_string_free (buffer, FALSE);
}


/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
