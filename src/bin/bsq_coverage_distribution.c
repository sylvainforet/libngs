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

static int  intcmp          (const int         *i1,
                             const int         *i2);

static void  rle_encode     (CallbackData      *data);

static void cleanup_data    (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData data;

  parse_args (&data, &argc, &argv);

  load_fasta (&data);
  load_bsq (&data);
  rle_encode (&data);

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
  qsort (data->reference_coverage,
         data->reference_size,
         sizeof (*data->reference_coverage),
         (int(*)(const void*, const void*))intcmp);
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

static int
intcmp (const int *i1,
        const int *i2)
{
  return *i1 - *i2;
}


static void
rle_encode (CallbackData *data)
{
  unsigned long int i;
  int               current;
  int               nb;

  current = -1;
  nb      = 0;
  for (i = 0; i < data->reference_size; i++)
    {
      if (data->reference_coverage[i] != current)
        {
          current = data->reference_coverage[i];
          if (nb > 0)
            {
              g_print ("%d\n", nb);
              nb = 0;
            }
          g_print ("%d\t", current);
        }
      nb++;
    }
  if (nb > 0)
    g_print ("%d\n", nb);
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
