/**
 *
 */

#include <stdlib.h>
#include <string.h>

#include "ngs_bsq.h"
#include "ngs_fasta.h"

/*********************/
/* ReferenceHashData */
/*********************/

typedef struct _ReferenceHashData ReferenceHashData;

struct _ReferenceHashData
{
  unsigned long int size;
  unsigned long int offset;
};

ReferenceHashData* reference_hash_data_new  (void);

void               reference_hash_data_free (ReferenceHashData *data);

/****************/
/* CallbackData */
/****************/

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  char             *input_fasta;
  char             *input_bsq;

  GHashTable       *reference_hash;
  GString          *reference_seq;
  int              *reference_coverage;

  unsigned long int reference_size;
};

static void parse_args      (CallbackData      *data,
                             int               *argc,
                             char            ***argv);

static void load_fasta      (CallbackData      *data);

static int  iter_fasta_func (FastaSeq          *seq,
                             CallbackData      *data);

static void load_bsq        (CallbackData      *data);

static int  iter_bsq_func   (BsqRecord         *rec,
                             CallbackData      *data);

static void get_content     (CallbackData      *data);

static void cleanup_data    (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData data;

  parse_args (&data, &argc, &argv);

  load_fasta (&data);
  load_bsq (&data);
  get_content (&data);
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
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  context = g_option_context_new ("FASTA  BSQ - Distribution of coverage");
  g_option_context_add_group (context, get_bsq_option_group ());
  g_option_context_add_main_entries (context, entries, NULL);
  if (!g_option_context_parse (context, argc, argv, &error))
    {
      g_printerr ("[ERROR] Option parsing failed: %s\n", error->message);
      exit (1);
    }
  g_option_context_free (context);

  if (*argc < 2)
    {
      g_printerr ("[ERROR] No input fasta file provided\n");
      exit (1);
    }
  if (*argc < 3)
    {
      g_printerr ("[ERROR] No input bsq file provided\n");
      exit (1);
    }
  data->input_fasta        = (*argv)[1];
  data->input_bsq          = (*argv)[2];
  data->reference_hash     = NULL;
  data->reference_coverage = NULL;
}

static void
load_fasta (CallbackData *data)
{
  GError *error = NULL;
  size_t  size;

  data->reference_size = 0;
  data->reference_hash = g_hash_table_new_full (g_str_hash,
                                                g_str_equal,
                                                g_free,
                                                (GDestroyNotify)reference_hash_data_free);
  data->reference_seq  = g_string_sized_new (1024 * 1024);
  iter_fasta (data->input_fasta,
              (FastaIterFunc)iter_fasta_func,
              data,
              &error);

  if (error)
    {
      g_printerr ("[ERROR] Loading fasta sequences failed: %s\n", error->message);
      g_error_free (error);
      cleanup_data (data);
      exit (1);
    }

  size                     = data->reference_size * sizeof (*data->reference_coverage);
  data->reference_coverage = g_malloc (size);
  memset (data->reference_coverage, 0, size);
}


static int
iter_fasta_func (FastaSeq     *seq,
                 CallbackData *data)
{
  ReferenceHashData *hash_data;

  hash_data             = reference_hash_data_new ();
  hash_data->size       = seq->size;
  hash_data->offset     = data->reference_size;
  data->reference_size += seq->size;
  /* TODO make sure the key is actually unique */
  g_hash_table_insert (data->reference_hash,
                       strdup (seq->name),
                       hash_data);
  data->reference_seq   = g_string_append (data->reference_seq,
                                           seq->seq);
  return 1;
}

static void
load_bsq (CallbackData *data)
{
  GError *error = NULL;

  iter_bsq (data->input_bsq,
            (BsqIterFunc)iter_bsq_func,
            data,
            &error);
  if (error)
    {
      g_printerr ("[ERROR] Iterating bsq failed: %s\n", error->message);
      g_error_free (error);
      cleanup_data (data);
      exit (1);
    }
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
      ReferenceHashData *hash_data;
      int               *start;
      int                i;

      hash_data = g_hash_table_lookup (data->reference_hash, rec->ref);
      if (!hash_data)
        {
          g_printerr ("[WARNING] Reference %s not found\n", rec->ref);
          return 1;
        }

      start = data->reference_coverage + hash_data->offset + rec->loc - 1;
      for (i = 0; i < rec->size; i++)
        start[i]++;
    }
  return 1;
}

static void
get_content (CallbackData *data)
{
  unsigned long int n_covered_A       = 0;
  unsigned long int n_covered_T       = 0;
  unsigned long int n_covered_G       = 0;
  unsigned long int n_covered_C       = 0;
  unsigned long int n_covered_N       = 0;
  unsigned long int n_covered_tot     = 0;
  unsigned long int n_non_covered_A   = 0;
  unsigned long int n_non_covered_T   = 0;
  unsigned long int n_non_covered_G   = 0;
  unsigned long int n_non_covered_C   = 0;
  unsigned long int n_non_covered_N   = 0;
  unsigned long int n_non_covered_tot = 0;
  unsigned long int i                 = 0;

  for (i = 0; i < data->reference_size; i++)
    {
      const char c = data->reference_seq->str[i];

      if (data->reference_coverage[i])
        {
          switch (c)
            {
              case 'A':
                  n_covered_A++;
                  break;
              case 'T':
                  n_covered_T++;
                  break;
              case 'G':
                  n_covered_G++;
                  break;
              case 'C':
                  n_covered_C++;
                  break;
              case 'N':
                  n_covered_N++;
                  break;
              default:
                  g_printerr ("[WARNING] Unknown character at position %ld: %c\n", i, c);
                  break;
            }
        }
      else
        {
          switch (c)
            {
              case 'A':
                  n_non_covered_A++;
                  break;
              case 'T':
                  n_non_covered_T++;
                  break;
              case 'G':
                  n_non_covered_G++;
                  break;
              case 'C':
                  n_non_covered_C++;
                  break;
              case 'N':
                  n_non_covered_N++;
                  break;
              default:
                  g_printerr ("[WARNING] Unknown character at position %ld: %c\n", i, c);
                  break;
            }
        }
    }
  n_covered_tot     = n_covered_A +
                      n_covered_T +
                      n_covered_G +
                      n_covered_C +
                      n_covered_N;
  n_non_covered_tot = n_non_covered_A +
                      n_non_covered_T +
                      n_non_covered_G +
                      n_non_covered_C +
                      n_non_covered_N;
  g_print ("Composition of covered bases:\n"
           " A: %9ld (%6.2f)\n"
           " T: %9ld (%6.2f)\n"
           " G: %9ld (%6.2f)\n"
           " C: %9ld (%6.2f)\n"
           " N: %9ld (%6.2f)\n",
           n_covered_A, 100. * n_covered_A / n_covered_tot,
           n_covered_T, 100. * n_covered_T / n_covered_tot,
           n_covered_G, 100. * n_covered_G / n_covered_tot,
           n_covered_C, 100. * n_covered_C / n_covered_tot,
           n_covered_N, 100. * n_covered_N / n_covered_tot);
  g_print ("Composition of non covered bases:\n"
           " A: %9ld (%6.2f)\n"
           " T: %9ld (%6.2f)\n"
           " G: %9ld (%6.2f)\n"
           " C: %9ld (%6.2f)\n"
           " N: %9ld (%6.2f)\n",
           n_non_covered_A, 100. * n_non_covered_A / n_non_covered_tot,
           n_non_covered_T, 100. * n_non_covered_T / n_non_covered_tot,
           n_non_covered_G, 100. * n_non_covered_G / n_non_covered_tot,
           n_non_covered_C, 100. * n_non_covered_C / n_non_covered_tot,
           n_non_covered_N, 100. * n_non_covered_N / n_non_covered_tot);
}

static void
cleanup_data (CallbackData *data)
{
  if (data->reference_hash)
    {
      g_hash_table_destroy (data->reference_hash);
      data->reference_hash = NULL;
    }
  if (data->reference_coverage)
    {
      g_free (data->reference_coverage);
      data->reference_coverage = NULL;
    }
  if (data->reference_seq)
    {
      g_string_free (data->reference_seq, TRUE);
      data->reference_seq = NULL;
    }
}

/*********************/
/* ReferenceHashData */
/*********************/

ReferenceHashData*
reference_hash_data_new (void)
{
  ReferenceHashData *data;

  data = g_slice_new0 (ReferenceHashData);

  return data;
}

void
reference_hash_data_free (ReferenceHashData *data)
{

  if (data)
    g_slice_free (ReferenceHashData, data);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
