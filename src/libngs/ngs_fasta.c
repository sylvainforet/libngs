/**
 *
 */

#include <stdlib.h>
#include <string.h>

#include "fasta.h"
#include "fasta_flex.h"


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
iter_fasta (char         *path,
            FastaIterFunc func,
            void         *data,
            GError      **error)
{
  if (fasta_parser_name == NULL || g_strcmp0 (fasta_parser_name, "flex") == 0)
    iter_fasta_flex (path, func, data, error);
  else
    {
      g_printerr ("[ERROR] Unknown fasta parser: %s\n",
                  fasta_parser_name);
      /* TODO raise an error instead of exiting */
      exit (1);
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

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
