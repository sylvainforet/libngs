/**
 *
 */

#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "ngs_binseq.h"
#include "ngs_fasta.h"
#include "ngs_fastq.h"
#include "ngs_kmerhash.h"
#include "ngs_memalloc.h"
#include "ngs_utils.h"

/* fast itoa implementarion, does not zero terminate the buffer and only works
 * with positive numbers */
#define uitoa_no0(i, buf, buf_size, ret, n_chars)   \
do {                                                \
  ret = buf + buf_size;                             \
  do {                                              \
      *--ret = '0' + (i % 10);                      \
      i /= 10;                                      \
      ++n_chars;                                    \
  } while (i != 0);                                 \
} while (0)

typedef struct _CallbackData CallbackData;

typedef int (*IterCharSeqFunc) (CallbackData        *data,
                                const char          *seq,
                                unsigned long int    size);

struct _CallbackData
{
  char              *input_path;
  char              *output_path;
  GIOChannel        *output_channel;

  IterCharSeqFunc    iter_char_seq_func;

  KmerHashTable     *htable;
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

static int               kmer_equal_default         (const unsigned char *key1,
                                                     const unsigned char *key2);

static unsigned long int kmer_hash_default          (const unsigned char *key);

/* Optimised Functions */

static int               kmer_equal_k16             (const unsigned char *key1,
                                                     const unsigned char *key2);

static unsigned long int kmer_hash_k16              (const unsigned char *key);

static int               iter_char_seq_k16          (CallbackData        *data,
                                                     const char          *seq,
                                                     unsigned long int    size);

static int               kmer_equal_k32             (const unsigned char *key1,
                                                     const unsigned char *key2);

static unsigned long int kmer_hash_k32              (const unsigned char *key);

static int               iter_char_seq_k32          (CallbackData        *data,
                                                     const char          *seq,
                                                     unsigned long int    size);

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
    kmer_hash_table_destroy (data.htable);
  if (data.output_path)
    g_free (data.output_path);
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

  context = g_option_context_new ("FILE - Count the number of kmers in a fasta/fastq file");
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
  data->input_path = (*argv)[1];

  word_size_bytes  = data->k_bytes;

  if (data->k == 16)
    {
      data->htable = kmer_hash_table_new ((GHashFunc)kmer_hash_k16,
                                          (GEqualFunc)kmer_equal_k16,
                                          data->k_bytes);
      data->iter_char_seq_func = iter_char_seq_k16;
    }
  else if (data->k > 16 && data->k < 32)
    {
      /* TODO this could be optimised a little further.  It seems that when k
       * gets close to 32, performance drops significantly */
      data->htable = kmer_hash_table_new ((GHashFunc)kmer_hash_k16,
                                          (GEqualFunc)kmer_equal_default,
                                          data->k_bytes);
      data->iter_char_seq_func = iter_char_seq_default;
    }
  else if (data->k == 32)
    {
      data->htable = kmer_hash_table_new ((GHashFunc)kmer_hash_k32,
                                          (GEqualFunc)kmer_equal_k32,
                                          data->k_bytes);
      data->iter_char_seq_func = iter_char_seq_k32;
    }
  else if (data->k > 32)
    {
      data->htable = kmer_hash_table_new ((GHashFunc)kmer_hash_k32,
                                          (GEqualFunc)kmer_equal_default,
                                          data->k_bytes);
      data->iter_char_seq_func = iter_char_seq_default;
    }
  else
    {
      data->htable = kmer_hash_table_new ((GHashFunc)kmer_hash_default,
                                          (GEqualFunc)kmer_equal_default,
                                          data->k_bytes);
      data->iter_char_seq_func = iter_char_seq_default;
    }

  if (data->verbose)
    g_printerr ("Parsing %s with k = %d\n", data->input_path, data->k);
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

  ret = data->iter_char_seq_func (data, seq, size);
  if (ret != 1)
    return ret;
  if (data->do_revcomp)
    ret = data->iter_char_seq_func (data, rev_comp_in_place (seq, size), size);

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
      kmer_hash_table_insert (data->htable, data->tmp_kmer);
    }
  return 1;
}

static void
print_results (CallbackData *data)
{
  KmerHashTableIter  iter;
  char              *buffer;
  char              *buffer_ptr;
  KmerHashNode      *hnode;
  GError            *error     = NULL;
  int               use_stdout = 1;

  if (data->verbose)
    g_printerr ("Printing results\n");
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

#define DIGITS_BUFF_SPACE 32
  buffer = g_alloca (data->k + DIGITS_BUFF_SPACE);
  buffer[data->k + DIGITS_BUFF_SPACE - 1] = '\n';
  kmer_hash_table_iter_init (&iter, data->htable);
  while ((hnode = kmer_hash_table_iter_next (&iter)) != NULL)
    {
      int j;
      if (data->bin_out)
        {
          for (j = 0; j < data->k_bytes; j++)
            buffer[j] = hnode->kmer[j];
          *((unsigned long int*)(buffer + j)) = GULONG_TO_BE (hnode->count);
          g_io_channel_write_chars (data->output_channel,
                                    (char*)buffer, data->k_bytes + sizeof (unsigned long int),
                                    NULL, &error);
        }
      else
        {
          unsigned long int n;
          j             = 1;
          n             = hnode->count;
          uitoa_no0 (n, buffer, DIGITS_BUFF_SPACE - 1 + data->k, buffer_ptr, j);
          *--buffer_ptr = ' ';
          j            += 1 + data->k;
          buffer_ptr   -= data->k;
          bin_to_char_prealloc (buffer_ptr, hnode->kmer, data->k);
          g_io_channel_write_chars (data->output_channel,
                                    buffer_ptr, j,
                                    NULL, &error);
        }
      if (error)
        {
          g_printerr ("[ERROR] Writing to output file `%s' failed: %s\n",
                      data->output_path,
                      error->message);
          exit (1);
        }
    }
#undef DIGITS_BUFF_SPACE

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

static int
iter_char_seq_k32 (CallbackData        *data,
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
      if (k == 32)
        {
          j++;
          break;
        }
    }
  if (k < 32)
    return 1;
  char_to_bin_prealloc (data->tmp_kmer, seq + j - k, data->k);
  data->tmp_kmer[0] <<= BITS_PER_NUC;
  data->tmp_kmer[1] <<= BITS_PER_NUC;
  data->tmp_kmer[2] <<= BITS_PER_NUC;
  data->tmp_kmer[3] <<= BITS_PER_NUC;
  data->tmp_kmer[4] <<= BITS_PER_NUC;
  data->tmp_kmer[5] <<= BITS_PER_NUC;
  data->tmp_kmer[6] <<= BITS_PER_NUC;
  data->tmp_kmer[7] <<= BITS_PER_NUC;
  for (i = j - k; i < maxi; i++)
    {
      char c;

      /* Skip over Ns */
      c = seq[i + data->k - 1];
      if (c == 'N' || c == 'n')
        {
          k = 0;
          for (j = i + 32; j < size; j++)
            {
              if (seq[j] == 'N' || seq[j] == 'n')
                k = 0;
              else
                k++;
              if (k == 32)
                {
                  j++;
                  break;
                }
            }
          if (k < 32)
            return 1;
          i = j - k;
          char_to_bin_prealloc (data->tmp_kmer, seq + i, data->k);
        }
      else
        {
          data->tmp_kmer[0] >>= BITS_PER_NUC;
          data->tmp_kmer[1] >>= BITS_PER_NUC;
          data->tmp_kmer[2] >>= BITS_PER_NUC;
          data->tmp_kmer[3] >>= BITS_PER_NUC;
          data->tmp_kmer[4] >>= BITS_PER_NUC;
          data->tmp_kmer[5] >>= BITS_PER_NUC;
          data->tmp_kmer[6] >>= BITS_PER_NUC;
          data->tmp_kmer[7] >>= BITS_PER_NUC;

          data->tmp_kmer[0] |= char_to_bin_table[(unsigned char)seq[i +  3]] << (BITS_PER_BYTE - BITS_PER_NUC);
          data->tmp_kmer[1] |= char_to_bin_table[(unsigned char)seq[i +  7]] << (BITS_PER_BYTE - BITS_PER_NUC);
          data->tmp_kmer[2] |= char_to_bin_table[(unsigned char)seq[i + 11]] << (BITS_PER_BYTE - BITS_PER_NUC);
          data->tmp_kmer[3] |= char_to_bin_table[(unsigned char)seq[i + 15]] << (BITS_PER_BYTE - BITS_PER_NUC);
          data->tmp_kmer[4] |= char_to_bin_table[(unsigned char)seq[i + 19]] << (BITS_PER_BYTE - BITS_PER_NUC);
          data->tmp_kmer[5] |= char_to_bin_table[(unsigned char)seq[i + 23]] << (BITS_PER_BYTE - BITS_PER_NUC);
          data->tmp_kmer[6] |= char_to_bin_table[(unsigned char)seq[i + 27]] << (BITS_PER_BYTE - BITS_PER_NUC);
          data->tmp_kmer[7] |= char_to_bin_table[(unsigned char)seq[i + 31]] << (BITS_PER_BYTE - BITS_PER_NUC);
        }
      kmer_hash_table_insert (data->htable, data->tmp_kmer);
    }
  return 1;
}

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
      if (k == 16)
        {
          j++;
          break;
        }
    }
  if (k < 16)
    return 1;
  char_to_bin_prealloc (data->tmp_kmer, seq + j - k, data->k);
  data->tmp_kmer[0] <<= BITS_PER_NUC;
  data->tmp_kmer[1] <<= BITS_PER_NUC;
  data->tmp_kmer[2] <<= BITS_PER_NUC;
  data->tmp_kmer[3] <<= BITS_PER_NUC;
  for (i = j - k; i < maxi; i++)
    {
      char c;

      /* Skip over Ns */
      c = seq[i + data->k - 1];
      if (c == 'N' || c == 'n')
        {
          k = 0;
          for (j = i + 16; j < size; j++)
            {
              if (seq[j] == 'N' || seq[j] == 'n')
                k = 0;
              else
                k++;
              if (k == 16)
                {
                  j++;
                  break;
                }
            }
          if (k < 16)
            return 1;
          i = j - k;
          char_to_bin_prealloc (data->tmp_kmer, seq + i, data->k);
        }
      else
        {
          data->tmp_kmer[0] >>= BITS_PER_NUC;
          data->tmp_kmer[1] >>= BITS_PER_NUC;
          data->tmp_kmer[2] >>= BITS_PER_NUC;
          data->tmp_kmer[3] >>= BITS_PER_NUC;

          data->tmp_kmer[0] |= char_to_bin_table[(unsigned char)seq[i +  3]] << (BITS_PER_BYTE - BITS_PER_NUC);
          data->tmp_kmer[1] |= char_to_bin_table[(unsigned char)seq[i +  7]] << (BITS_PER_BYTE - BITS_PER_NUC);
          data->tmp_kmer[2] |= char_to_bin_table[(unsigned char)seq[i + 11]] << (BITS_PER_BYTE - BITS_PER_NUC);
          data->tmp_kmer[3] |= char_to_bin_table[(unsigned char)seq[i + 15]] << (BITS_PER_BYTE - BITS_PER_NUC);
        }
      kmer_hash_table_insert (data->htable, data->tmp_kmer);
    }
  return 1;
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

static int
kmer_equal_k16 (const unsigned char *key1,
                const unsigned char *key2)
{
  const guint32 i1 = *(guint32*)key1;
  const guint32 i2 = *(guint32*)key2;
  return i1 == i2;
}

static unsigned long int
kmer_hash_k16 (const unsigned char *key)
{
  return (*(guint32*)key) * 2654435769U;;
}

static int
kmer_equal_k32 (const unsigned char *key1,
                const unsigned char *key2)
{
  const guint64 i1 = *(guint64*)key1;
  const guint64 i2 = *(guint64*)key2;
  return i1 == i2;
}

static unsigned long int
kmer_hash_k32 (const unsigned char *key)
{
  guint32 *k = (guint32*) key;
  return (k[0] ^ k[1]) * 2654435769U;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
