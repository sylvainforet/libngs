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
  unsigned long int a_qual[N_QUAL];
  unsigned long int t_qual[N_QUAL];
  unsigned long int g_qual[N_QUAL];
  unsigned long int c_qual[N_QUAL];
  unsigned long int n_qual[N_QUAL];
};

static int  iter_func             (FastqSeq        *fastq,
                                   LetterQualData  *data);

static void init_letter_qual_data (LetterQualData  *data);

static void print_letter_qual     (LetterQualData  *data);

int
main (int    argc,
      char **argv)
{
  LetterQualData  data;
  GError         *error = NULL;

  if (argc < 2)
    {
      g_printerr ("Usage: letter_qual FILE\n");
      return 1;
    }

  init_letter_qual_data (&data);
  iter_fastq (argv[1],
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

static int
iter_func  (FastqSeq       *fastq,
            LetterQualData *data)
{
  unsigned int i;

  for (i = 0; fastq->seq[i] && fastq->qual[i]; i++)
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
