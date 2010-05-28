/**
 *
 */

#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "ngs_binseq.h"
#include "ngs_fasta.h"

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char          *input_path;

  GHashTable    *htable;
  GArray        *keys;
  GArray        *values;
  unsigned char *tmp_key;

  unsigned int   bin_kmer_bytes;
  int            k;
};

/* For use in the general-purpose comparison function */
static size_t bin_kmer_bytes;

static void parse_args (CallbackData      *data,
                        int               *argc,
                        char            ***argv);

static int  iter_func  (FastaSeq          *fasta,
                        CallbackData      *data);

static unsigned int bin_seq_hash_8   (const void *key);
static unsigned int bin_seq_hash_16  (const void *key);
static unsigned int bin_seq_hash_32  (const void *key);
static unsigned int bin_seq_hash_64  (const void *key);
static unsigned int bin_seq_hash_xx  (const void *key);

static int          bin_seq_equal_8  (const void *a,
                                      const void *b);
static int          bin_seq_equal_16 (const void *a,
                                      const void *b);
static int          bin_seq_equal_32 (const void *a,
                                      const void *b);
static int          bin_seq_equal_64 (const void *a,
                                      const void *b);
static int          bin_seq_equal_xx (const void *a,
                                      const void *b);

/*
static int          compare_int      (const void *a,
                                      const void *b);
                                      */

int
main (int    argc,
      char **argv)
{
  CallbackData  data;
  GError       *error = NULL;
  unsigned int  i;

  parse_args (&data, &argc, &argv);

  iter_fasta (data.input_path,
              (FastaIterFunc)iter_func,
              &data,
              &error);

  if (error)
    {
      g_printerr ("[ERROR] Iterating sequences failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
      return 1;
    }
  for (i = 0; i < data.values->len; i++)
    {
      char *tmp = bin_to_char ((unsigned char*)data.keys->data + data.bin_kmer_bytes * i, data.k);
      g_print ("%s %d\n", tmp, g_array_index (data.values, int, i));
      g_free (tmp);
    }

  g_hash_table_destroy (data.htable);
  g_array_free (data.keys, TRUE);
  g_array_free (data.values, TRUE);
  g_free (data.tmp_key);

  return 0;
}

static void
parse_args (CallbackData      *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      {"kmer", 'k', 0, G_OPTION_ARG_INT, &data->k, "K-mer size", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;
  GHashFunc       hash_func;
  GEqualFunc      equal_func;

  data->k = 8;

  context = g_option_context_new ("FILE - Count the number of kmers in a fasta file");
  g_option_context_add_group (context, get_fasta_option_group ());
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
  if (data->k < 1)
    {
      g_printerr ("[ERROR] The kmer size must be larger than 0\n");
      exit (1);
    }
  else if (data->k <= 4)
    {
      hash_func  = bin_seq_hash_8;
      equal_func = bin_seq_equal_8;
    }
  else if (data->k <= 8)
    {
      hash_func  = bin_seq_hash_16;
      equal_func = bin_seq_equal_16;
    }
  else if (data->k <= 16)
    {
      hash_func  = bin_seq_hash_32;
      equal_func = bin_seq_equal_32;
    }
  else if (data->k <= 32)
    {
      hash_func  = bin_seq_hash_64;
      equal_func = bin_seq_equal_64;
    }
  else
    {
      hash_func  = bin_seq_hash_xx;
      equal_func = bin_seq_equal_xx;
    }
  data->input_path      = (*argv)[1];
  data->htable          = g_hash_table_new (hash_func, equal_func);
  data->bin_kmer_bytes  = (data->k + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;
  data->tmp_key         = g_malloc0 (data->bin_kmer_bytes);
  data->keys            = g_array_new (FALSE, FALSE, data->bin_kmer_bytes);
  data->values          = g_array_new (FALSE, FALSE, sizeof (int));
  bin_kmer_bytes        = data->bin_kmer_bytes;
}

static int
iter_func (FastaSeq     *fasta,
           CallbackData *data)
{
  const unsigned long int maxi = fasta->size - data->k + 1;
  unsigned long int       i;

  for (i = 0; i < maxi; i++)
    {
      gpointer     val = NULL;
      gpointer     orig_key = NULL;
      unsigned int j;

      /* Clear the key */
      for (j = 0; j < data->bin_kmer_bytes; j++)
        data->tmp_key[j] = 0;

      /* TODO handle sequences with Ns */
      char_to_bin_prealloc (data->tmp_key, fasta->seq + i, data->k);
      if (g_hash_table_lookup_extended (data->htable, data->tmp_key, &orig_key, &val))
        {
          /*
          for (j = 0; j < data->values->len; j++)
            g_print ("%d ", g_array_index (data->values, int, j));
          g_print ("\n");

          g_print (">>> %ld %s %d %p\n", i, bin_to_char (orig_key, data->k), *((int*)val), val);
          */
          g_array_index (data->values, int, GPOINTER_TO_INT (val))++;
        }
      else
        {
          int new_val = 1;

          data->keys   = g_array_append_vals (data->keys, data->tmp_key, 1);
          data->values = g_array_append_vals (data->values, &new_val, 1);
          g_hash_table_insert (data->htable,
                               data->keys->data   + (data->keys->len   - 1) * data->bin_kmer_bytes,
                               GINT_TO_POINTER (data->values->len - 1));

          /*
          for (j = 0; j < data->values->len; j++)
            g_print ("%d ", g_array_index (data->values, int, j));
          g_print ("\n");

          g_print (">>> --- %ld %s %d\n", i,
                   bin_to_char ((unsigned char*)data->keys->data + (data->keys->len - 1) * data->bin_kmer_bytes, data->k),
                   *(int*)(data->values->data + (data->values->len - 1) * sizeof (int)));
          */
        }
    }
  return 1;
}

static unsigned int
bin_seq_hash_8 (const void *key)
{
  return *(guint8*)key;
}

static unsigned int
bin_seq_hash_16 (const void *key)
{
  return *(guint16*)key;
}

static unsigned int
bin_seq_hash_32 (const void *key)
{
  return *(guint32*)key;
}

static unsigned int
bin_seq_hash_64 (const void *key)
{
  return *(guint32*)key;
}

static unsigned int
bin_seq_hash_xx (const void *key)
{
  return *(guint32*)key;
}

static int
bin_seq_equal_8 (const void *a,
                 const void *b)
{
  return ((*(guint8*)a) - (*(guint8*)b)) == 0;
}

static int
bin_seq_equal_16 (const void *a,
                  const void *b)
{
  return ((*(guint16*)a) - (*(guint16*)b)) == 0;
}

static int
bin_seq_equal_32 (const void *a,
                  const void *b)
{
  return ((*(guint32*)a) - (*(guint32*)b)) == 0;
}

static int
bin_seq_equal_64 (const void *a,
                  const void *b)
{
  return ((*(guint64*)a) - (*(guint64*)b)) == 0;
}

static int
bin_seq_equal_xx (const void *a,
                  const void *b)
{
  return memcmp (a, b, bin_kmer_bytes) == 0;
}

/*
static int
compare_int (const void *a,
             const void *b)
{
  return (*(int*)a) - (*(int*)b);
}
*/

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
