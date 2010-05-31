/**
 *
 */

#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "ngs_binseq.h"
#include "ngs_fasta.h"
#include "ngs_fastq.h"
#include "ngs_utils.h"


typedef struct _KmerHashTable KmerHashTable;

typedef unsigned long int (*KmerHashTableHashFunc)  (KmerHashTable *htable,
                                                     unsigned char *key);

typedef int               (*KmerHashTableEqualFunc) (KmerHashTable *htable,
                                                     unsigned char *key1,
                                                     unsigned char *key2);

struct _KmerHashTable
{
  KmerHashTableHashFunc  hash_func;
  KmerHashTableEqualFunc equal_func;

  GPtrArray       **buckets;
  GStringChunk     *kwords;

  unsigned long int size;
  unsigned long int nb_elems;
  unsigned int      word_size_letters;
  unsigned int      word_size_bytes;
};

typedef struct _KmerHashNode KmerHashNode;

struct _KmerHashNode
{
  unsigned char    *kmer;
  unsigned long int n;
};

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char              *input_path;
  char              *output_path;
  GIOChannel        *output_channel;

  KmerHashTable     *htable;
  char              *tmp_str;

  unsigned long int *harray;

  unsigned long int  hsize;
  unsigned long int  n_seqs;
  unsigned long int  freq_report;
  unsigned int       k;
  int                do_revcomp;
  int                is_fastq;
  int                bin_out;
  int                verbose;
};

static KmerHashTable*    kmer_hash_table_new        (unsigned long int    size,
                                                     unsigned int         word_size_letters);

static void              kmer_hash_table_free       (KmerHashTable       *htable);

/*
static KmerHashNode*     kmer_hash_table_lookup     (KmerHashTable       *htable,
                                                     unsigned char       *key);

static void              kmer_hash_table_insert     (KmerHashTable       *htable,
                                                     unsigned char       *key);
*/

static void              kmer_hash_table_inc        (KmerHashTable       *htable,
                                                     unsigned char       *key);

static int               kmer_equal_default         (KmerHashTable       *htable,
                                                     unsigned char       *key1,
                                                     unsigned char       *key2);

static unsigned long int kmer_hash_default          (KmerHashTable       *htable,
                                                     unsigned char       *key);

static KmerHashNode*     kmer_hash_node_new         (void);

static void              kmer_hash_node_free        (KmerHashNode        *hnode);

static void              parse_args                 (CallbackData        *data,
                                                     int                 *argc,
                                                     char              ***argv);

static int               iter_func_fasta            (FastaSeq            *fasta,
                                                     CallbackData        *data);

static int               iter_func_fastq            (FastqSeq            *fasta,
                                                     CallbackData        *data);

static int               iter_func_char             (char                *seq,
                                                     unsigned long int    size,
                                                     CallbackData        *data);

static int               iter_char_seq_kx           (CallbackData        *data,
                                                     const char          *seq,
                                                     unsigned long int    size);

static void              kmer_hash_node_print       (KmerHashNode        *hnode,
                                                     CallbackData        *data);

static void              print_results              (CallbackData        *data);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;
  GError       *error = NULL;

  parse_args (&data, &argc, &argv);

  if (data.is_fastq)
    iter_fastq (data.input_path,
                (FastqIterFunc)iter_func_fastq,
                &data,
                &error);
  else
    iter_fasta (data.input_path,
                (FastaIterFunc)iter_func_fasta,
                &data,
                &error);

  print_results (&data);

  if (error)
    {
      g_printerr ("[ERROR] Iterating sequences failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
      return 1;
    }
  if (data.harray)
    g_free (data.harray);
  if (data.htable)
    kmer_hash_table_free (data.htable);
  if (data.tmp_str)
    g_free (data.tmp_str);
  if (data.output_path)
    g_free (data.output_path);

  return 0;
}

static void
parse_args (CallbackData      *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      {"out",      'o', 0, G_OPTION_ARG_FILENAME, &data->output_path, "Output file", NULL},
      {"hsize",    's', 0, G_OPTION_ARG_INT,      &data->hsize,       "Number of buckets in the hash table", NULL},
      {"kmer",     'k', 0, G_OPTION_ARG_INT,      &data->k,           "K-mer size", NULL},
      {"freqrep",  'e', 0, G_OPTION_ARG_INT,      &data->freq_report, "Verbose-report frequency", NULL},
      {"revcomp",  'r', 0, G_OPTION_ARG_NONE,     &data->do_revcomp,  "Also scan the reverse complement", NULL},
      {"fastq",    'q', 0, G_OPTION_ARG_NONE,     &data->is_fastq,    "Input is in fastq format", NULL},
      {"binout",   'b', 0, G_OPTION_ARG_NONE,     &data->bin_out,     "Write output in binary format", NULL},
      {"verbose",  'v', 0, G_OPTION_ARG_NONE,     &data->verbose,     "Verbose output", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->hsize       = 1024 * 1024;
  data->k           = 0;
  data->do_revcomp  = 0;
  data->is_fastq    = 0;
  data->bin_out     = 0;
  data->verbose     = 0;
  data->n_seqs      = 0;
  data->freq_report = 1000000;
  data->output_path = strdup("-");
  data->harray      = NULL;
  data->htable      = NULL;
  data->tmp_str     = NULL;

  context = g_option_context_new ("FILE - Count the number of kmers in a fasta file");
  g_option_context_add_group (context, get_fasta_option_group ());
  g_option_context_add_group (context, get_fastq_option_group ());
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

  if (data->hsize < 1)
    {
      g_printerr ("[ERROR] The number of buckets must be larger than 0\n");
      exit (1);
    }
  if (data->k < 1)
    {
      g_printerr ("[ERROR] The kmer size must be larger than 0\n");
      exit (1);
    }

  data->htable     = kmer_hash_table_new (data->hsize, data->k);
  if (data->verbose)
    g_printerr ("Created hashing table with %ld buckets and %d letters per word\n", data->hsize, data->k);
  data->input_path = (*argv)[1];
}

static int
iter_func_fasta (FastaSeq     *fasta,
                 CallbackData *data)
{
  return iter_func_char (fasta->seq, fasta->size, data);
}

static int
iter_func_fastq (FastqSeq     *fastq,
                 CallbackData *data)
{
  return iter_func_char (fastq->seq, fastq->size, data);
}

static int
iter_func_char (char              *seq,
                unsigned long int  size,
                CallbackData      *data)
{
  int            ret;

  ret = iter_char_seq_kx (data, seq, size);
  if (ret != 1)
    return ret;
  if (data->do_revcomp)
    ret = iter_char_seq_kx (data, rev_comp_in_place (seq, size), size);

  data->n_seqs++;
  if (data->verbose && data->n_seqs % data->freq_report == 0)
    g_printerr ("Parsed %ld sequences\n", data->n_seqs);

  return ret;
}

static int
iter_char_seq_kx (CallbackData     *data,
                  const char       *seq,
                  unsigned long int size)
{
  unsigned long int       i;
  unsigned long int       j;
  unsigned long int       k;
  const unsigned long int maxi = size - data->k + 1;

  if (size < data->k)
    return 1;

  /* Find the first index of a kmer without Ns */
  for (j = 0, k = 0; j < size; j++)
    {
      if (seq[j] == 'N' || seq[j] == 'n')
        k = 0;
      else
        k++;
      if (k == data->k)
        {
          j++;
          break;
        }
    }
  if (k < data->k)
    return 1;

  for (i = j - k; i < maxi; i++)
    {
      char c;

      /* Skip over Ns */
      c = seq[i + data->k - 1];
      if (c == 'N' || c == 'n')
        {
          k = 0;
          for (j = i + data->k; j < size; j++)
            {
              if (seq[j] == 'N' || seq[j] == 'n')
                k = 0;
              else
                k++;
              if (k == data->k)
                {
                  j++;
                  break;
                }
            }
          if (k < data->k)
            return 1;
          i = j - k;
        }
      kmer_hash_table_inc (data->htable, (unsigned char*)seq + i);
    }
  return 1;
}

static void
print_results (CallbackData *data)
{
  GError            *error      = NULL;
  int               use_stdout = 1;
  unsigned long int i;
  unsigned long int j;

  /* Open */
  if (!data->output_path || !*data->output_path || (data->output_path[0] == '-' && data->output_path[1] == '\0'))
    data->output_channel = g_io_channel_unix_new (STDOUT_FILENO);
  else
    {
      use_stdout     = 0;
      data->output_channel = g_io_channel_new_file (data->output_path, "w", &error);
      if (error)
        {
          g_printerr ("[ERROR] failed to open output file `%s': %s\n",
                      data->output_path,
                      error->message);
          exit (1);
        }
    }
  g_io_channel_set_encoding (data->output_channel, NULL, NULL);

  for (i = 0; i < data->htable->size; i++)
    if (data->htable->buckets[i])
      for (j = 0; j < data->htable->buckets[i]->len; j++)
        kmer_hash_node_print ((KmerHashNode*)data->htable->buckets[i]->pdata[j], data);

  /* Close */
  if (!use_stdout)
    {
      g_io_channel_shutdown (data->output_channel, TRUE, &error);
      if (error)
        {
          g_printerr ("[ERROR] Closing output file `%s' failed: %s\n",
                      data->output_path,
                      error->message);
          g_error_free (error);
        }
    }
  g_io_channel_unref (data->output_channel);
}

static void
kmer_hash_node_print (KmerHashNode *hnode,
                      CallbackData *data)
{
  GError  *error = NULL;
  GString *buffer;

  buffer = g_string_new (NULL);
  if (data->bin_out)
    {
      const unsigned long int val = GULONG_TO_BE (hnode->n);
      g_string_append_len (buffer, (char*)hnode->kmer, data->k);
      g_string_append_len (buffer, (char*)&val, sizeof (unsigned long int));
      g_io_channel_write_chars (data->output_channel,
                                buffer->str,
                                data->k + sizeof (unsigned long int),
                                NULL,
                                &error);
    }
  else
    {
      g_string_append_len (buffer, (char*)hnode->kmer, data->k);
      g_string_append_printf (buffer, " %ld\n", hnode->n);
      g_io_channel_write_chars (data->output_channel,
                                buffer->str,
                                -1,
                                NULL,
                                &error);
    }
  g_string_free (buffer, TRUE);

  if (error)
    {
      g_printerr ("[ERROR] Writing to output file `%s' failed: %s\n",
                  data->output_path,
                  error->message);
      g_error_free (error);
      error = NULL;
    }
}

static KmerHashTable*
kmer_hash_table_new (unsigned long int size,
                     unsigned int      word_size_letters)
{
  KmerHashTable *htable;

  htable                    = g_slice_new0 (KmerHashTable);
  htable->size              = size;
  htable->word_size_letters = word_size_letters;
  htable->word_size_bytes   = htable->word_size_letters;
  htable->buckets           = g_malloc0 (htable->size * sizeof (*htable->buckets));
  htable->kwords            = g_string_chunk_new (word_size_letters);
  htable->hash_func         = kmer_hash_default;
  htable->equal_func        = kmer_equal_default;

  return htable;
}

static void
kmer_hash_table_free (KmerHashTable *htable)
{
  if (htable)
    {
      if (htable->kwords)
        g_string_chunk_free (htable->kwords);
      if (htable->buckets)
        {
          unsigned long int i;

          for (i = 0; i < htable->size; i++)
            if (htable->buckets[i])
              g_ptr_array_free (htable->buckets[i], TRUE);

          g_free (htable->buckets);
        }
      g_slice_free (KmerHashTable, htable);
    }
}

static unsigned long int
kmer_hash_default (KmerHashTable *htable,
                   unsigned char *key)
{
  unsigned long int hash = 5381;
  unsigned int      i    = 0;

  do
    hash = (hash * 33) ^ key[i];
  while (++i < htable->word_size_bytes);

  return hash % htable->size;
}

static int
kmer_equal_default (KmerHashTable *htable,
                    unsigned char *key1,
                    unsigned char *key2)
{
  unsigned int i;

  for (i = 0; i < htable->word_size_bytes; i++)
    if (key1[i] != key2[i])
      return 0;

  return 1;
}

/*
static KmerHashNode*
kmer_hash_table_lookup (KmerHashTable *htable,
                        unsigned char *key)
{
  GPtrArray        *bucket;
  unsigned long int idx;
  unsigned long int i;

  idx    = htable->hash_func (htable, key);
  bucket = htable->buckets[idx];
  if (!bucket)
    return NULL;
  for (i = 0; i < bucket->len; i++)
    if (htable->equal_func (htable, key, ((KmerHashNode*)bucket->pdata[i])->kmer))
      return (KmerHashNode*)bucket->pdata[i];

  return NULL;
}

static void
kmer_hash_table_insert (KmerHashTable *htable,
                        unsigned char *key)
{
  KmerHashNode     *hnode;
  GPtrArray        *bucket;
  unsigned long int idx;
  unsigned long int i;

  idx    = htable->hash_func (htable, key);
  bucket = htable->buckets[idx];
  if (bucket)
    for (i = 0; i < bucket->len; i++)
      if (htable->equal_func (htable, key, ((KmerHashNode*)bucket->pdata[i])->kmer))
        return;
  htable->buckets[idx] = g_ptr_array_new_with_free_func ((GDestroyNotify)kmer_hash_node_free);
  hnode                = kmer_hash_node_new ();
  hnode->kmer          = (unsigned char*)g_string_chunk_insert_len (htable->kwords, (char*)key, htable->word_size_letters);
  hnode->n             = 0;
  g_ptr_array_add (htable->buckets[idx], hnode);
  htable->nb_elems++;
}
*/

static void
kmer_hash_table_inc (KmerHashTable *htable,
                     unsigned char *key)
{
  KmerHashNode     *hnode;
  GPtrArray        *bucket;
  unsigned long int idx;
  unsigned long int i;

  idx    = htable->hash_func (htable, key);
  bucket = htable->buckets[idx];
  if (bucket)
    {
    for (i = 0; i < bucket->len; i++)
      if (htable->equal_func (htable, key, ((KmerHashNode*)bucket->pdata[i])->kmer))
        {
          ((KmerHashNode*)bucket->pdata[i])->n++;
          return;
        }
    }
  else
    htable->buckets[idx] = g_ptr_array_new_with_free_func ((GDestroyNotify)kmer_hash_node_free);

  hnode       = kmer_hash_node_new ();
  hnode->kmer = (unsigned char*)g_string_chunk_insert_len (htable->kwords, (char*)key, htable->word_size_letters);
  hnode->n    = 1;
  g_ptr_array_add (htable->buckets[idx], hnode);
  htable->nb_elems++;
}

static KmerHashNode*
kmer_hash_node_new (void)
{
  KmerHashNode *hnode;

  hnode = g_slice_new0 (KmerHashNode);

  return hnode;
}

static void
kmer_hash_node_free (KmerHashNode *hnode)
{
  if (hnode)
    g_slice_free (KmerHashNode, hnode);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
