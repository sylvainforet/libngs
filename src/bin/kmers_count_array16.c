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

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char              *input_path;
  char              *output_path;
  GIOChannel        *output_channel;

  unsigned long int *array;
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

static int               iter_char_seq_k16          (CallbackData        *data,
                                                     const char          *seq,
                                                     unsigned long int    size);

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
      /* {"kmer",     'k', 0, G_OPTION_ARG_INT,      &data->k,           "K-mer size", NULL}, */
      {"freqrep",  'e', 0, G_OPTION_ARG_INT,      &data->freq_report, "Verbose-report frequency", NULL},
      {"revcomp",  'r', 0, G_OPTION_ARG_NONE,     &data->do_revcomp,  "Also scan the reverse complement", NULL},
      {"fastq",    'q', 0, G_OPTION_ARG_NONE,     &data->is_fastq,    "Input is in fastq format", NULL},
      {"binout",   'b', 0, G_OPTION_ARG_NONE,     &data->bin_out,     "Write output in binary format", NULL},
      {"verbose",  'v', 0, G_OPTION_ARG_NONE,     &data->verbose,     "Verbose output", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->k           = 16;
  data->do_revcomp  = 0;
  data->is_fastq    = 0;
  data->bin_out     = 0;
  data->verbose     = 0;
  data->n_seqs      = 0;
  data->freq_report = 1000000;
  data->output_path = strdup("-");
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

  data->k_bytes    = (data->k + NUCS_PER_BYTE - 1) / NUCS_PER_BYTE;
  data->tmp_kmer   = g_malloc0 (data->k_bytes);
  data->input_path = (*argv)[1];
  word_size_bytes  = data->k_bytes;

  if (data->verbose)
    g_printerr ("Allocating array of size %lu\n", (1L<<32));
  data->array      = g_malloc0 ((1L<<32) * sizeof (*data->array));
  if (data->verbose)
    g_printerr ("DONE !!!\n");

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

  ret = iter_char_seq_k16 (data, seq, size);
  if (ret != 1)
    return ret;
  if (data->do_revcomp)
    ret = iter_char_seq_k16 (data, rev_comp_in_place (seq, size), size);

  data->n_seqs++;
  if (data->verbose && data->n_seqs % data->freq_report == 0)
    g_printerr ("Parsed %ld sequences\n", data->n_seqs);

  return ret;
}

static inline void
guint32_to_char  (char                *dest,
                  guint32              src)
{
  dest[ 0]   = bin_to_char_table[src       & 3];
  dest[ 1]   = bin_to_char_table[src >>  2 & 3];
  dest[ 2]   = bin_to_char_table[src >>  4 & 3];
  dest[ 3]   = bin_to_char_table[src >>  6 & 3];
  dest[ 4]   = bin_to_char_table[src >>  8 & 3];
  dest[ 5]   = bin_to_char_table[src >> 10 & 3];
  dest[ 6]   = bin_to_char_table[src >> 12 & 3];
  dest[ 7]   = bin_to_char_table[src >> 14 & 3];
  dest[ 8]   = bin_to_char_table[src >> 16 & 3];
  dest[ 9]   = bin_to_char_table[src >> 18 & 3];
  dest[10]   = bin_to_char_table[src >> 20 & 3];
  dest[11]   = bin_to_char_table[src >> 22 & 3];
  dest[12]   = bin_to_char_table[src >> 24 & 3];
  dest[13]   = bin_to_char_table[src >> 26 & 3];
  dest[14]   = bin_to_char_table[src >> 28 & 3];
  dest[15]   = bin_to_char_table[src >> 30 & 3];
}

static void
print_results (CallbackData *data)
{
  struct {
      guint32 i;
      unsigned long int n;
  } bin_out;
  unsigned long int  i;
  char              *buffer;
  char              *tmp_str;
  char              *itoa_ptr;
  GError            *error      = NULL;
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

  tmp_str    = g_alloca (data->k);
  buffer     = g_alloca (64);
  buffer[63] = '\n';
  for (i = 0; i < (1L<<32); i++)
    {
      if (!data->array[i])
        continue;
      if (data->bin_out)
        {
          /* FIXME whatchout for endianness */
          bin_out.i = (guint32)i;
          bin_out.n = data->array[i];
          g_io_channel_write_chars (data->output_channel,
                                    (char*)&bin_out, sizeof (bin_out),
                                    NULL, NULL);
        }
      else
        {
          unsigned long int n;
          int               j = 1;

          guint32_to_char (tmp_str, (guint32)i);
          g_io_channel_write_chars (data->output_channel,
                                    tmp_str, data->k,
                                    NULL, NULL);
          n = data->array[i];
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

static inline  guint32
char16_to_guint32 (const char *src)
{
  guint32 ret = 0;

  ret  |= char_to_bin_table[(unsigned char)src[15]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[14]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[13]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[12]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[11]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[10]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[ 9]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[ 8]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[ 7]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[ 6]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[ 5]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[ 4]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[ 3]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[ 2]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[ 1]];
  ret <<= 2;
  ret  |= char_to_bin_table[(unsigned char)src[ 0]];
  return ret;
}

static int
iter_char_seq_k16 (CallbackData        *data,
                   const char          *seq,
                   unsigned long int    size)
{
  unsigned long int       idx;
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
  idx   = char16_to_guint32 (seq + j - k);
  idx <<= 2;
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
          idx = char16_to_guint32 (seq + j - k);
        }
      else
        {
          idx >>= 2;
          idx  |= char_to_bin_table[(unsigned char)seq[i + 15]] << 30;
        }
      data->array[idx]++;
    }
  return 1;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
