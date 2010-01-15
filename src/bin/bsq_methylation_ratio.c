/**
 *
 */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ngs_bsq.h"
#include "ngs_fasta.h"
#include "ngs_fastq.h"
#include "ngs_seq_db.h"
#include "ngs_utils.h"


typedef struct _MethData MethData;

struct _MethData
{
  unsigned int n_meth;
  unsigned int n_unmeth;

  /*
  unsigned int n_quals
  char        *quals;
  */
};

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char              *ref_path;
  char             **fastq_paths;
  char             **bsq_paths;
  char              *output_path;
  char              *add_path;

  SeqDB             *reads;
  SeqDB             *ref;
  MethData         **meth_index;
  MethData          *meth_data;

  int                min_qual;
  int                verbose;
  int                print_letter;
  int                print_all;
};

static void parse_args    (CallbackData      *data,
                           int               *argc,
                           char            ***argv);

static void load_data     (CallbackData      *data);

static void map_data      (CallbackData      *data);

static int  iter_bsq_func (BsqRecord         *rec,
                           CallbackData      *data);

static void cleanup_data  (CallbackData      *data);

static void print_meth    (CallbackData      *data);

static void load_previous (CallbackData      *data);


int
main (int    argc,
      char **argv)
{
  CallbackData  data;

  parse_args (&data, &argc, &argv);
  if (data.verbose)
    g_print (">>> Loading Data\n");
  load_data (&data);
  if (data.add_path)
    load_previous (&data);
  if (data.verbose)
    g_print (">>> Mapping Data\n");
  map_data (&data);
  if (data.verbose)
    g_print (">>> Writing Data\n");
  print_meth (&data);
  cleanup_data (&data);

  return 0;
}

static void
parse_args (CallbackData      *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      {"reference", 'r', 0, G_OPTION_ARG_FILENAME,       &data->ref_path,     "Reference genome file", NULL},
      {"bsq",       'b', 0, G_OPTION_ARG_FILENAME_ARRAY, &data->bsq_paths,    "Bsmap bsq file(s)", NULL},
      {"fastq",     'f', 0, G_OPTION_ARG_FILENAME_ARRAY, &data->fastq_paths,  "Fastq file(s)", NULL},
      {"out",       'o', 0, G_OPTION_ARG_FILENAME,       &data->output_path,  "Output file", NULL},
      {"add",       'a', 0, G_OPTION_ARG_FILENAME,       &data->add_path,     "Add results to this file", NULL},
      {"min_qual",  'm', 0, G_OPTION_ARG_INT,            &data->min_qual,     "Minimum base quality", NULL},
      {"verbose",   'v', 0, G_OPTION_ARG_NONE,           &data->verbose,      "Verbose output", NULL},
      {"letter",    'l', 0, G_OPTION_ARG_NONE,           &data->print_letter, "Prepend a column with the letter", NULL},
      {"all",       'w', 0, G_OPTION_ARG_NONE,           &data->print_all,    "Prints all positions (implies l)", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->ref_path     = NULL;
  data->fastq_paths  = NULL;
  data->bsq_paths    = NULL;
  data->output_path  = strdup("-");
  data->add_path     = NULL;
  data->min_qual     = 0;
  data->verbose      = 0;
  data->print_letter = 0;
  data->print_all    = 0;

  context = g_option_context_new (" - Finds all the potentially methylated Cs");
  g_option_context_add_group (context, get_fasta_option_group ());
  g_option_context_add_group (context, get_fastq_option_group ());
  g_option_context_add_group (context, get_bsq_option_group ());
  g_option_context_add_main_entries (context, entries, NULL);
  if (!g_option_context_parse (context, argc, argv, &error))
    {
      g_printerr ("[ERROR] Option parsing failed: %s\n", error->message);
      exit (1);
    }
  g_option_context_free (context);

  if (!data->ref_path)
    {
      g_printerr ("[ERROR] You must specify a reference genome with -r\n");
      exit (1);
    }
  if (!data->fastq_paths || !*data->fastq_paths)
    {
      g_printerr ("[ERROR] You must specify one or more fastq files with -f\n");
      exit (1);
    }
  if (!data->bsq_paths || !*data->fastq_paths)
    {
      g_printerr ("[ERROR] You must specify one or more bsq files with -b\n");
      exit (1);
    }
  if (data->min_qual >= 35)
    g_printerr ("[WARNING] Minimum quality is set to %d, "
                "it is likely that no base will be considered\n",
                data->min_qual);
  if (data->print_all)
    data->print_letter = 1;
}

static void
load_data (CallbackData *data)
{
  GError         *error = NULL;
  char          **tmp;
  unsigned long   n_cg;
  unsigned long   i;

  data->ref   = seq_db_new ();
  data->reads = seq_db_new ();

  if (data->verbose)
    g_print (">>> Loading Fasta\n");
  seq_db_load_fasta (data->ref, data->ref_path, &error);
  if (error)
    {
      g_printerr ("[ERROR] Loading reference failed: %s\n", error->message);
      exit (1);
    }
  if (data->verbose)
    g_print (">>> Loading Bsq\n");
  for (tmp = data->fastq_paths; *tmp; tmp++)
    {
      seq_db_load_fastq (data->reads, *tmp, &error);
      if (error)
        {
          g_printerr ("[ERROR] Loading reads failed: %s\n", error->message);
          exit (1);
        }
    }
  n_cg = 0;
  for (i = 0; i < data->ref->total_size; i++)
    if (data->ref->seqs[i] == 'C' || data->ref->seqs[i] == 'G')
      ++n_cg;

  data->meth_index = g_malloc0 (data->ref->total_size * sizeof (*data->meth_index) );
  data->meth_data  = g_malloc0 (n_cg * sizeof (*data->meth_data) );

  n_cg = 0;
  for (i = 0; i < data->ref->total_size; i++)
    if (data->ref->seqs[i] == 'C' || data->ref->seqs[i] == 'G')
      data->meth_index[i] = &data->meth_data[n_cg++];
}

static int
iter_bsq_func (BsqRecord    *rec,
               CallbackData *data)
{
  int name_len;

  name_len = strlen (rec->name);
  if (rec->flag == BSQ_MAP_UM &&
      ((rec->name[name_len - 2] == '/' &&
       rec->name[name_len - 1] == '1' &&
       (rec->strand == BSQ_STRAND_W ||
        rec->strand == BSQ_STRAND_C)) ||
      (rec->name[name_len - 2] == '/' &&
       rec->name[name_len - 1] == '2' &&
       (rec->strand == BSQ_STRAND_WC ||
        rec->strand == BSQ_STRAND_CC))))
    {
      SeqDBElement *read_elem;
      SeqDBElement *ref_elem;
      char         *read;
      char         *qual;
      char         *ref;
      unsigned int  start_ref;
      unsigned int  i;
      int           is_ref_rev = 0;

      read_elem = g_hash_table_lookup (data->reads->index, rec->name);
      if (!read_elem)
        {
          g_printerr ("[WARNING] Read `%s' not found\n", rec->name);
          return 1;
        }
      ref_elem = g_hash_table_lookup (data->ref->index, rec->ref);
      if (!ref_elem)
        {
          g_printerr ("[WARNING] Reference `%s' not found\n", rec->ref);
          return 1;
        }

      start_ref = ref_elem->offset + rec->loc - 1;
      read      = data->reads->seqs + read_elem->offset;
      qual      = data->reads->quals + read_elem->offset;
      ref       = data->ref->seqs + start_ref;

      switch (rec->strand)
        {
          case BSQ_STRAND_W:
              /* Dont reverse */
              break;
          case BSQ_STRAND_C:
              /* Reverse ref */
              ref = rev_comp_in_place (ref, read_elem->size);
              is_ref_rev = 1;
              break;
          case BSQ_STRAND_WC:
              /* Reverse read */
              read = rev_comp_in_place (read, read_elem->size);
              qual = rev_in_place (qual, read_elem->size);
              break;
          case BSQ_STRAND_CC:
              /* Reverse both read and ref */
              ref  = rev_comp_in_place (ref, read_elem->size);
              qual = rev_in_place (qual, read_elem->size);
              read = rev_comp_in_place (read, read_elem->size);
              is_ref_rev = 1;
              break;
          default:
              break;
        }

      for (i = 0; i < read_elem->size; i++)
        {
          if (ref[i] == 'C')
            {
              if (qual[i] >= data->min_qual)
                {
                  unsigned int meth_idx;

                  meth_idx = i;
                  if (is_ref_rev)
                    meth_idx = read_elem->size - i - 1;
                  if (read[i] == 'C')
                    data->meth_index[start_ref + meth_idx]->n_meth++;
                  else if (read[i] == 'T')
                    data->meth_index[start_ref + meth_idx]->n_unmeth++;
                }
            }
        }

      /* Put the reference back into place
       * Should not need to have the reads back in the right orientation ...
       */
      if (is_ref_rev)
        ref  = rev_comp_in_place (ref, read_elem->size);
    }
  return 1;
}

static void
map_data (CallbackData *data)
{
  char  **tmp;

  for (tmp = data->bsq_paths; *tmp; tmp++)
    {
      GError *error = NULL;

      iter_bsq (*tmp,
                (BsqIterFunc)iter_bsq_func,
                data,
                &error);
      if (error)
        {
          g_printerr ("[ERROR] Loading bsq file `%s' failed: %s\n",
                      *tmp, error->message);
          exit (1);
        }
    }
}

static void
cleanup_data (CallbackData *data)
{
  if (!data)
    return;
  if (data->fastq_paths)
    g_strfreev (data->fastq_paths);
  if (data->bsq_paths)
    g_strfreev (data->bsq_paths);
  if (data->ref_path)
    g_free (data->ref_path);
  if (data->output_path)
    g_free (data->output_path);
  if (data->meth_index)
    g_free (data->meth_index);
  if (data->meth_data)
    g_free (data->meth_data);
  seq_db_free (data->ref);
  seq_db_free (data->reads);
}

static void
print_meth (CallbackData *data)
{
  GHashTableIter iter;
  GIOChannel    *channel;
  SeqDBElement  *elem;
  GError        *error      = NULL;
  int            use_stdout = 1;

  /* Open */
  if (data->output_path[0] == '-' && data->output_path[1] == '\0')
    channel = g_io_channel_unix_new (STDOUT_FILENO);
  else
    {
      use_stdout = 0;
      channel = g_io_channel_new_file (data->output_path, "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] Opening output file failed: %s\n", error->message);
          exit (1);
        }
    }

  g_hash_table_iter_init (&iter, data->ref->index);
  while (g_hash_table_iter_next (&iter, NULL, (gpointer*)&elem))
    {
      GString      *buffer;
      MethData    **ref_meth;
      unsigned long i;

      ref_meth = data->meth_index + elem->offset;
      buffer   = g_string_new (NULL);
      g_string_printf (buffer, ">%s\n", elem->name);
      for (i = 0; i < elem->size; i++)
        {
        if (ref_meth[i])
          {
            if (data->print_letter)
              g_string_append_printf (buffer,
                                      "%c\t%ld\t%u\t%u\n",
                                      data->ref->seqs[elem->offset + i],
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
        else if (data->print_all)
          g_string_append_printf (buffer,
                                  "%c\n",
                                  data->ref->seqs[elem->offset + i]);
        }
      g_io_channel_write_chars (channel,
                                buffer->str,
                                buffer->len,
                                NULL,
                                &error);
      if (error)
        {
          g_printerr ("[ERROR] Writing to output file failed: %s\n",
                      error->message);
          g_error_free (error);
          error = NULL;
        }

      g_string_free (buffer, TRUE);
    }

  /* Close */
  if (!use_stdout)
    {
      g_io_channel_shutdown (channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing output file failed: %s\n",
                      error->message);
          g_error_free (error);
        }
    }
  g_io_channel_unref (channel);
}

static void
load_previous (CallbackData *data)
{
  GIOChannel    *channel;
  SeqDBElement  *elem  = NULL;
  GError        *error = NULL;
  char          *line  = NULL;
  gsize          length;
  gsize          endl;

  /* Open */
  channel = g_io_channel_new_file (data->add_path, "r", &error);
  if (error)
    {
      g_printerr ("[ERROR] Opening previous results file `%s' failed: %s\n",
                  data->add_path,
                  error->message);
      exit (1);
    }

  /* Parse */
  while (G_IO_STATUS_NORMAL == g_io_channel_read_line (channel, &line, &length, &endl, &error))
    {
      line[endl] = '\0';
      if (line[0] == '>')
        {
          elem = g_hash_table_lookup (data->ref->index, line + 1);
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
              data->meth_index[elem->offset + offset]->n_meth   = g_ascii_strtoll (values[1], NULL, 10);
              data->meth_index[elem->offset + offset]->n_unmeth = g_ascii_strtoll (values[2], NULL, 10);
            }
          else if (n_fields == 4)
            {
              offset = g_ascii_strtoll (values[1], NULL, 10);
              data->meth_index[elem->offset + offset]->n_meth   = g_ascii_strtoll (values[2], NULL, 10);
              data->meth_index[elem->offset + offset]->n_unmeth = g_ascii_strtoll (values[3], NULL, 10);
            }
          else if (n_fields != 1)
            g_printerr ("[WARNING] Could not parse previous run's line\n");

          g_strfreev (values);
        }
      g_free (line);
      line = NULL;
    }
  if (error)
    {
      g_printerr ("[ERROR] Reading previous results file failed: %s\n",
                  error->message);
      exit (1);
    }
  if (line)
    g_free (line);

  /* Close */
  g_io_channel_shutdown (channel, TRUE, &error);
  if (error)
    {
      g_printerr ("[ERROR] Closing previous results file failed: %s\n",
                  error->message);
      g_error_free (error);
    }
  g_io_channel_unref (channel);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
