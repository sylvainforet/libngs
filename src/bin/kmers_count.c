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

/* fast itoa implementarion, does not zero terminate the buffer and only works
 * with strictly positive numbers */
#define uitoa_no0(i, buf, buf_size, ret, n_chars)   \
do {                                                \
  ret = buf + buf_size;                             \
  do {                                              \
      *--ret = '0' + (i % 10);                      \
      i /= 10;                                      \
      ++n_chars;                                    \
  } while (i != 0);                                 \
} while (0)

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

  GHashTable        *htable;
  GStringChunk      *chunks;
  unsigned char     *tmp_kmer;

  unsigned long int  n_seqs;
  unsigned long int  freq_report;
  unsigned int       k;
  unsigned int       k_bytes;
  int                do_revcomp;
  int                is_fastq;
  int                bin_out;
  int                verbose;
};

static unsigned int word_size_bytes;

/* Default Functions */

static int               kmer_equal_default         (const unsigned char *key1,
                                                     const unsigned char *key2);

static unsigned long int kmer_hash_default          (const unsigned char *key);

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

static int               iter_char_seq_default      (CallbackData        *data,
                                                     const char          *seq,
                                                     unsigned long int    size);

static void              print_results              (CallbackData        *data);

/* Optimised Functions */

static int               kmer_equal_k16             (const unsigned char *key1,
                                                     const unsigned char *key2);

static unsigned long int kmer_hash_k16              (const unsigned char *key);

#if 0
static int               iter_char_seq_k16          (CallbackData        *data,
                                                     const char          *seq,
                                                     unsigned long int    size);
#endif

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
  if (data.htable)
    g_hash_table_destroy (data.htable);
  if (data.output_path)
    g_free (data.output_path);
  if (data.chunks)
    g_string_chunk_free (data.chunks);
  if (data.tmp_kmer)
    g_free (data.tmp_kmer);

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

  data->k           = 0;
  data->do_revcomp  = 0;
  data->is_fastq    = 0;
  data->bin_out     = 0;
  data->verbose     = 0;
  data->n_seqs      = 0;
  data->freq_report = 1000000;
  data->output_path = strdup("-");
  data->htable      = NULL;
  data->tmp_kmer    = NULL;

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

  if (data->k < 1)
    {
      g_printerr ("[ERROR] The kmer size must be larger than 0\n");
      exit (1);
    }

  data->k_bytes    = (data->k + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;
  data->tmp_kmer   = g_malloc0 (data->k_bytes);
  data->chunks     = g_string_chunk_new (data->k_bytes);
  data->input_path = (*argv)[1];

  word_size_bytes  = data->k_bytes;

  if (data->k == 16)
    data->htable = g_hash_table_new_full ((GHashFunc)kmer_hash_k16,
                                          (GEqualFunc)kmer_equal_k16,
                                          NULL,
                                          (GDestroyNotify)kmer_hash_node_free);
  else
    data->htable = g_hash_table_new_full ((GHashFunc)kmer_hash_default,
                                          (GEqualFunc)kmer_equal_default,
                                          NULL,
                                          (GDestroyNotify)kmer_hash_node_free);
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

  ret = iter_char_seq_default (data, seq, size);
  if (ret != 1)
    return ret;
  if (data->do_revcomp)
    ret = iter_char_seq_default (data, rev_comp_in_place (seq, size), size);

  data->n_seqs++;
  if (data->verbose && data->n_seqs % data->freq_report == 0)
    g_printerr ("Parsed %ld sequences\n", data->n_seqs);

  return ret;
}

static int
iter_char_seq_default (CallbackData     *data,
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
      KmerHashNode *node;

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
      char_to_bin_prealloc (data->tmp_kmer, seq + i, data->k);
      if (g_hash_table_lookup_extended (data->htable, (char*)data->tmp_kmer, NULL, (void**)&node))
        node->n++;
      else
        {
          node       = kmer_hash_node_new ();
          node->n    = 1;
          node->kmer = (unsigned char*)g_string_chunk_insert_len (data->chunks, (char*)data->tmp_kmer, data->k_bytes);
          g_hash_table_insert (data->htable, node->kmer, node);
        }
    }
  return 1;
}

static void
print_results (CallbackData *data)
{
  GHashTableIter     iter;
  unsigned char     *key;
  char              *buffer;
  char              *tmp_str;
  char              *itoa_ptr;
  KmerHashNode      *hnode;
  GError            *error      = NULL;
  int               use_stdout = 1;

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
  g_io_channel_set_buffer_size (data->output_channel, 1024 * 1024 * 32);

  tmp_str    = g_alloca (data->k);
  buffer     = g_alloca (64);
  buffer[63] = '\n';
  g_hash_table_iter_init (&iter, data->htable);
  while (g_hash_table_iter_next (&iter, (void**)&key, (void**)&hnode))
    {
      if (data->bin_out)
        {
          const unsigned long int val = GULONG_TO_BE (hnode->n);
          g_io_channel_write_chars (data->output_channel,
                                    (char*)hnode->kmer, data->k_bytes,
                                    NULL, NULL);
          g_io_channel_write_chars (data->output_channel,
                                    (char*)&val, sizeof (unsigned long int),
                                    NULL, NULL);
        }
      else
        {
          unsigned long int n;
          int               j = 1;

          bin_to_char_prealloc (tmp_str, hnode->kmer, data->k);
          g_io_channel_write_chars (data->output_channel,
                                    tmp_str, data->k,
                                    NULL, NULL);
          n = hnode->n;
          uitoa_no0 (n, buffer, 63, itoa_ptr, j);
          *--itoa_ptr = ' ';
          ++j;
          g_io_channel_write_chars (data->output_channel,
                                    itoa_ptr, j,
                                    NULL, NULL);
        }
    }

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

static unsigned long int
kmer_hash_default (const unsigned char *key)
{
  /*
  unsigned long int hash = 5381;
  unsigned int      i    = 0;

  do
    hash = (hash * 33) ^ key[i];
  while (++i < word_size_bytes);
  */
  /*
  unsigned long int hash = 0;
  unsigned int      i    = 0;

  do
    hash = (hash << 5) - hash + key[i];
  while (++i < word_size_bytes);
  */
  unsigned long int hash = 0;
  unsigned int      i    = 0;
  do
    {
      hash <<= 8;
      hash  |= key[i];
    }
  while (++i < word_size_bytes);

  return ((unsigned int)(hash * 2654435769U));;
}

static int
kmer_equal_default (const unsigned char *key1,
                    const unsigned char *key2)
{
  unsigned int i;

  for (i = 0; i < word_size_bytes; i++)
    if (key1[i] != key2[i])
      return 0;

  return 1;
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

/* Optimised Functions */

static int
kmer_equal_k16 (const unsigned char *key1,
                const unsigned char *key2)
{
  int b1 = (key1[0] == key2[0]);
  int b2 = (key1[1] == key2[1]);
  int b3 = (key1[2] == key2[2]);
  int b4 = (key1[3] == key2[3]);
  b1 = b1 && b2;
  b3 = b3 && b4;
  return b1 && b3;
}

static unsigned long int
kmer_hash_k16 (const unsigned char *key)
{
  unsigned long int hash = key[3];
  hash <<= 8;
  hash  |= key[2];
  hash <<= 8;
  hash  |= key[1];
  hash <<= 8;
  hash  |= key[0];
  return ((unsigned int)(hash * 2654435769U));;
}

#if 0
static int
iter_char_seq_k16 (CallbackData        *data,
                   const char          *seq,
                   unsigned long int    size)
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
      KmerHashNode *node;

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
      char_to_bin_prealloc (data->tmp_kmer, seq + i, data->k);
      if (g_hash_table_lookup_extended (data->htable, (char*)data->tmp_kmer, NULL, (void**)&node))
        node->n++;
      else
        {
          node       = kmer_hash_node_new ();
          node->n    = 1;
          node->kmer = (unsigned char*)g_string_chunk_insert_len (data->chunks, (char*)data->tmp_kmer, data->k_bytes);
          g_hash_table_insert (data->htable, node->kmer, node);
        }
    }
  return 1;
}
#endif

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
