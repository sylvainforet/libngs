/**
 *
 */

#include <stdlib.h>
#include <unistd.h>

#include "fastq.h"

#define N_QUAL 40

typedef struct _LetterQualData LetterQualData;

struct _LetterQualData
{
  unsigned long int *quals[128];

  unsigned long int  a_qual[N_QUAL];
  unsigned long int  t_qual[N_QUAL];
  unsigned long int  g_qual[N_QUAL];
  unsigned long int  c_qual[N_QUAL];
  unsigned long int  n_qual[N_QUAL];

  char              *input_path;

  unsigned int       fast;
};

static void parse_args            (LetterQualData    *data,
                                   int               *argc,
                                   char            ***argv);

static int  iter_func             (FastqSeq          *fastq,
                                   LetterQualData    *data);

static int  iter_func_fast        (FastqSeq          *fastq,
                                   LetterQualData    *data);

static void init_letter_qual_data (LetterQualData    *data);

static void print_letter_qual     (LetterQualData    *data);

int
main (int    argc,
      char **argv)
{
  LetterQualData  data;
  GError         *error = NULL;

  parse_args (&data, &argc, &argv);

  if (data.fast)
    iter_fastq (data.input_path,
                (FastqIterFunc)iter_func_fast,
                &data,
                &error);
  else
    iter_fastq (data.input_path,
                (FastqIterFunc)iter_func,
                &data,
                &error);

  if (error)
    {
      g_printerr ("[ERROR] Iterating sequences failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
      return 1;
    }
  print_letter_qual (&data);

  return 0;
}

static void
parse_args (LetterQualData    *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      {"fast", 'f', 0, G_OPTION_ARG_NONE, &data->fast, "Faster and less robust method (use at your own risk)", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->fast = 0;

  context = g_option_context_new ("FILE - Letter-wise distribution of quality calls");
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
  data->input_path = (*argv)[1];

  init_letter_qual_data (data);
}

/**
 * Slower but more robust
 */

static int
iter_func (FastqSeq       *fastq,
           LetterQualData *data)
{
  unsigned int i;

  for (i = 0; i < fastq->size; i++)
    {
      unsigned long int *table = NULL;

      switch (fastq->seq[i])
        {
          case 'A':
              table = data->a_qual;
              break;
          case 'T':
              table = data->t_qual;
              break;
          case 'G':
              table = data->g_qual;
              break;
          case 'C':
              table = data->c_qual;
              break;
          case 'N':
              table = data->n_qual;
              break;
          default:
              g_printerr ("[WARNING] Unknown sequence character %c\n",
                          fastq->seq[i]);
              break;
        }
      if (table)
        {
          const unsigned int idx = fastq->qual[i] - 64;
          if (idx >= N_QUAL)
            g_printerr ("[WARNING] Unknown quality character %c\n",
                        fastq->qual[i]);
          else
            table[idx]++;
        }
    }
  return 1;
}

/**
 * Faster but no less robust
 */

static int
iter_func_fast (FastqSeq       *fastq,
                LetterQualData *data)
{
  unsigned int i;

  for (i = 0; i < fastq->size; i++)
    {
      const unsigned int letter = fastq->seq[i];
      const unsigned int qual   = fastq->qual[i] - 64;
      data->quals[letter][qual]++;
    }

  return 1;
}

static void
init_letter_qual_data (LetterQualData  *data)
{
  int i;

  for (i = 0; i < N_QUAL; i++)
    data->a_qual[i] = 0;
  for (i = 0; i < N_QUAL; i++)
    data->t_qual[i] = 0;
  for (i = 0; i < N_QUAL; i++)
    data->g_qual[i] = 0;
  for (i = 0; i < N_QUAL; i++)
    data->c_qual[i] = 0;
  for (i = 0; i < N_QUAL; i++)
    data->n_qual[i] = 0;

  if (data->fast)
    {
      data->quals['A'] = data->a_qual;
      data->quals['T'] = data->t_qual;
      data->quals['G'] = data->g_qual;
      data->quals['C'] = data->c_qual;
      data->quals['N'] = data->n_qual;
    }
}

static void
print_letter_qual (LetterQualData *data)
{
  int i;

  /* Header */
  g_print ("A\t"
           "T\t"
           "G\t"
           "C\t"
           "N\n");
  for (i = 0; i < 40; i++)
    {
      g_print ("%ld\t"
               "%ld\t"
               "%ld\t"
               "%ld\t"
               "%ld\n",
               data->a_qual[i],
               data->t_qual[i],
               data->g_qual[i],
               data->c_qual[i],
               data->n_qual[i]);
    }
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
