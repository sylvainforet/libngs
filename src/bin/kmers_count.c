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


typedef struct _CallbackData CallbackData;

typedef int (*IterBinSeqFunc) (CallbackData     *data,
                               unsigned char    *bin,
                               unsigned long int size);

static unsigned long int* make_hash_table_k8  (void);

static unsigned long int* make_hash_table_k16 (void);

static int                iter_bin_seq_k8     (CallbackData     *data,
                                               unsigned char    *bin,
                                               unsigned long int size);

static int                iter_bin_seq_k16    (CallbackData     *data,
                                               unsigned char    *bin,
                                               unsigned long int size);

static void               print_results_k8    (CallbackData     *data);

static void               print_results_k16   (CallbackData     *data);

typedef struct _BinSeqHashNode BinSeqHashNode;

struct _BinSeqHashNode
{
  unsigned char    *word;
  unsigned long int n;
};

static BinSeqHashNode* bin_seq_hash_node_new   (void);

static void            bin_seq_hash_node_free  (BinSeqHashNode      *hnode);

static void            bin_seq_hash_node_print (unsigned char       *key,
                                                BinSeqHashNode      *node,
                                                CallbackData        *data);

static unsigned int    bin_seq_hash            (const unsigned char *bin);

static int             bin_seq_equal           (const unsigned char *bin1,
                                                const unsigned char *bin2);

static int             iter_char_seq_kx        (CallbackData        *data,
                                                const char          *seq,
                                                unsigned long int    size);

static void            print_results_kx        (CallbackData        *data);

static unsigned int kmer_size;
static unsigned int kmer_bytes;

struct _CallbackData
{
  char              *input_path;
  char              *output_path;
  GIOChannel        *output_channel;

  IterBinSeqFunc     iter_bin_func;

  GHashTable        *htable;
  GStringChunk      *kwords;
  unsigned char     *tmp_word;
  char              *tmp_str;

  unsigned long int *harray;

  unsigned long int  n_seqs;
  unsigned long int  freq_report;
  unsigned int       k;
  unsigned int       k_bytes;
  int                do_revcomp;
  int                is_fastq;
  int                bin_out;
  int                fast;
  int                verbose;
};

static void parse_args       (CallbackData      *data,
                              int               *argc,
                              char            ***argv);

static int  iter_func_fasta  (FastaSeq          *fasta,
                              CallbackData      *data);

static int  iter_func_fastq  (FastqSeq          *fasta,
                              CallbackData      *data);

static int  iter_func_char   (char              *seq,
                              unsigned long int  size,
                              CallbackData      *data);

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

  print_results_kx (&data);

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
    g_hash_table_destroy (data.htable);
  if (data.kwords)
    g_string_chunk_free (data.kwords);
  if (data.tmp_word)
    g_free (data.tmp_word);
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
      {"kmer",     'k', 0, G_OPTION_ARG_INT,      &data->k,           "K-mer size", NULL},
      {"freqrep",  'e', 0, G_OPTION_ARG_INT,      &data->freq_report, "Verbose-report frequency", NULL},
      {"revcomp",  'r', 0, G_OPTION_ARG_NONE,     &data->do_revcomp,  "Also scan the reverse complement", NULL},
      {"fastq",    'q', 0, G_OPTION_ARG_NONE,     &data->is_fastq,    "Input is in fastq format", NULL},
      {"binout",   'b', 0, G_OPTION_ARG_NONE,     &data->bin_out,     "Write output in binary format", NULL},
      {"fast",     'f', 0, G_OPTION_ARG_NONE,     &data->fast,        "Use optimisations if available", NULL},
      {"verbose",  'v', 0, G_OPTION_ARG_NONE,     &data->verbose,     "Verbose output", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->k           = 0;
  data->do_revcomp  = 0;
  data->is_fastq    = 0;
  data->bin_out     = 0;
  data->fast        = 0;
  data->verbose     = 0;
  data->n_seqs      = 0;
  data->freq_report = 1000000;
  data->output_path = strdup("-");
  data->harray      = NULL;
  data->htable      = NULL;
  data->tmp_word    = NULL;
  data->tmp_str     = NULL;
  data->kwords      = NULL;

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

  data->k_bytes = (data->k + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;
  if (data->fast)
    {
      if (data->k == 8)
        {
          data->iter_bin_func = iter_bin_seq_k8;
          data->harray        = make_hash_table_k8 ();
        }
      else if (data->k == 16)
        {
          data->iter_bin_func = iter_bin_seq_k16;
          data->harray        = make_hash_table_k16 ();
        }
      else
        {
          g_printerr ("[ERROR] No optimised functions available for k=%d\n", data->k);
          exit (1);
        }
    }
  else
    {
      data->iter_bin_func    = NULL;
      data->htable           = g_hash_table_new_full ((GHashFunc)bin_seq_hash,
                                                      (GEqualFunc)bin_seq_equal,
                                                      NULL,
                                                      (GDestroyNotify)bin_seq_hash_node_free);
      data->kwords           = g_string_chunk_new (data->k_bytes);
      kmer_size              = data->k;
      kmer_bytes             = data->k_bytes;
      data->tmp_word         = g_malloc (data->k_bytes);
      data->tmp_str          = g_malloc (data->k + 1);
      data->tmp_str[data->k] = '\0';
    }
  data->input_path = (*argv)[1];
}

static int
iter_func_char (char              *seq,
                unsigned long int  size,
                CallbackData      *data)
{
  unsigned char *bin;
  int            ret;

  if (data->iter_bin_func)
    {
      bin = char_to_bin (seq, size);
      ret = data->iter_bin_func (data, bin, size);
      g_free (bin);
      if (ret != 1)
        return ret;
      if (data->do_revcomp)
        {
          bin = char_to_bin (rev_comp_in_place (seq, size), size);
          ret = data->iter_bin_func (data, bin, size);
          g_free (bin);
        }
    }
  else
    {
      ret = iter_char_seq_kx (data, seq, size);
      if (ret != 1)
        return ret;
      if (data->do_revcomp)
        ret = iter_char_seq_kx (data, rev_comp_in_place (seq, size), size);
    }
  data->n_seqs++;
  if (data->verbose && data->n_seqs % data->freq_report == 0)
    g_printerr ("Parsed %ld sequences\n", data->n_seqs);

  return ret;
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

static unsigned long int*
make_hash_table_k8 (void)
{
  unsigned long int *count_array;

  count_array = g_malloc0 ((1 << 16) * sizeof (*count_array));

  return count_array;
}

static unsigned long int*
make_hash_table_k16 (void)
{
  unsigned long int *count_array;

  count_array = g_malloc0 ((1L << 32) * sizeof (*count_array));

  return count_array;
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
      BinSeqHashNode *hnode;
      char            c;

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

      for (j = 0; j < data->k_bytes; j++)
        data->tmp_word[j] = 0;

      data->tmp_word = char_to_bin_prealloc (data->tmp_word, seq + i, data->k);
      if (g_hash_table_lookup_extended (data->htable, data->tmp_word, NULL, (gpointer*)&hnode))
        hnode->n++;
      else
        {
          hnode       = bin_seq_hash_node_new ();
          hnode->n    = 1;
          hnode->word = (unsigned char*)g_string_chunk_insert_len (data->kwords,
                                                                   (char*)data->tmp_word,
                                                                   data->k_bytes);
          g_hash_table_insert (data->htable, hnode->word, hnode);
        }
    }
  return 1;
}

/* FIXME this does not handle Ns */
static int
iter_bin_seq_k8 (CallbackData     *data,
                 unsigned char    *bin,
                 unsigned long int size)
{
  unsigned long int i;
  guint16           word;

  if (size < 8)
    return 1;
  
  word   = bin[1];
  word <<= 8;
  word  |= bin[0];
  data->harray[word]++;
  for (i = 8; i < size; i++)
    {
      unsigned long int address;
      unsigned int      offset;

      word  >>= BITS_PER_NUC;
      address = i / NUCS_PER_BYTE;
      offset  = (i % NUCS_PER_BYTE) * BITS_PER_NUC;
      word   |= ((bin[address] >> offset) & 3) << 14;
      data->harray[word]++;
    }

  return 1;
}

/* FIXME this does not handle Ns */
static int
iter_bin_seq_k16 (CallbackData     *data,
                  unsigned char    *bin,
                  unsigned long int size)
{
  unsigned long int i;
  guint32           word;

  if (size < 16)
    return 1;
  
  word   = bin[3];
  word <<= 8;
  word  |= bin[2];
  word <<= 8;
  word  |= bin[1];
  word <<= 8;
  word  |= bin[0];
  data->harray[word]++;
  for (i = 16; i < size; i++)
    {
      unsigned long int address;
      unsigned int      offset;

      word  >>= BITS_PER_NUC;
      address = i / NUCS_PER_BYTE;
      offset  = (i % NUCS_PER_BYTE) * BITS_PER_NUC;
      word   |= ((bin[address] >> offset) & 3) << 30;
      data->harray[word]++;
    }

  return 1;
}

static void
print_results_k8 (CallbackData *data)
{
  unsigned char  bin_buf[2];
  BinSeqHashNode hnode;
  unsigned int   i;

  for (i = 0; i < (1 << 16); i++)
    {
      if (!data->harray[i])
        continue;
      bin_buf[0] = (i & 255);
      bin_buf[1] = (i >> 8);

      hnode.word = bin_buf;
      hnode.n    = data->harray[i];
      bin_seq_hash_node_print (bin_buf, &hnode, data);
    }
}

static void
print_results_k16 (CallbackData *data)
{
  unsigned char  bin_buf[4];
  BinSeqHashNode hnode;
  unsigned int  i;

  for (i = 0; i < (1L << 32); i++)
    {
      if (!data->harray[i])
        continue;
      bin_buf[0] = (i & 255);
      bin_buf[1] = (i >> 8) & 255;
      bin_buf[2] = (i >> 16) & 255;
      bin_buf[3] = (i >> 24);

      hnode.word = bin_buf;
      hnode.n    = data->harray[i];
      bin_seq_hash_node_print (bin_buf, &hnode, data);
    }
}

static void
print_results_kx (CallbackData *data)
{
  GError *error      = NULL;
  int     use_stdout = 1;

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

  if (data->fast)
    {
      if (data->k == 8)
        print_results_k8 (data);
      else if (data->k == 16)
        print_results_k16 (data);
    }
  else
    g_hash_table_foreach (data->htable,
                          (GHFunc)bin_seq_hash_node_print,
                          data);
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

static BinSeqHashNode*
bin_seq_hash_node_new (void)
{
  BinSeqHashNode *hnode;

  hnode = g_slice_new0 (BinSeqHashNode);

  return hnode;
}

static void
bin_seq_hash_node_free (BinSeqHashNode *hnode)
{
  if (hnode)
    g_slice_free (BinSeqHashNode, hnode);
}

static unsigned int
bin_seq_hash (const unsigned char *bin)
{
  guint32      h = 0;
  unsigned int i;

  for (i = 0; i < kmer_bytes && i < 4; i++)
    {
      h <<= 8;
      h  |= bin[i];
      /*  This has one extra instruction, but keeps the first letters in the
       *  low bits
      h >>= 8;
      h |= bin[i] << 24;
      */
    }
  return h;
}

static int
bin_seq_equal (const unsigned char *bin1,
               const unsigned char *bin2)
{
  unsigned int i;

  for (i = 0 ; i < kmer_bytes; i++)
    if (bin1[i] != bin2[i])
      return 0;
  return 1;
}

static void
bin_seq_hash_node_print (unsigned char       *key,
                         BinSeqHashNode      *node,
                         CallbackData        *data)
{
  GError  *error = NULL;
  GString *buffer;

  buffer = g_string_new (NULL);
  if (data->bin_out)
    {
      const unsigned long int val = GULONG_TO_BE (node->n);
      g_string_append_len (buffer, (char*)node->word, data->k_bytes);
      g_string_append_len (buffer, (char*)&val, sizeof (unsigned long int));
      g_io_channel_write_chars (data->output_channel,
                                buffer->str,
                                data->k_bytes + sizeof (unsigned long int),
                                NULL,
                                &error);
    }
  else
    {
      data->tmp_str = bin_to_char_prealloc (data->tmp_str, node->word, data->k);
      g_string_printf (buffer, "%s %ld\n", data->tmp_str, node->n);
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

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
