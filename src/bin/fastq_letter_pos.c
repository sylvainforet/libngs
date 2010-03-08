/**
 *
 */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ngs_fastq.h"


#define N_LETTERS 5

#define NUC_A     0
#define NUC_T     1
#define NUC_G     2
#define NUC_C     3
#define NUC_N     4


typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  unsigned long int *counts[N_LETTERS];
  unsigned long int *ns;

  char              *input_path;

  unsigned int       current_size;
};

static void parse_args    (CallbackData      *data,
                           int               *argc,
                           char            ***argv);

static int  iter_func     (FastqSeq          *fastq,
                           CallbackData      *data);

static void print_letters (CallbackData      *data);

int
main (int    argc,
      char **argv)
{
  CallbackData  data;
  GError       *error = NULL;

  parse_args (&data, &argc, &argv);

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

  print_letters (&data);

  return 0;
}

static void
parse_args (CallbackData      *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      /*{"fast", 'f', 0, G_OPTION_ARG_NONE, &data->fast, "Faster and less robust method (use at your own risk)", NULL},*/
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  context = g_option_context_new ("FILE - Computes the letter composition at each position");
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
  data->input_path   = (*argv)[1];
  data->current_size = 0;
  data->ns           = NULL;
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  int i;

  if (data->current_size == 0)
    {
      data->ns = g_malloc (fastq->size * sizeof (*data->ns));
      data->ns = memset (data->ns, 0, fastq->size * sizeof (*data->ns));
      for (i = 0; i < N_LETTERS; i++)
        {
          data->counts[i] = g_malloc (fastq->size * sizeof (**data->counts));
          data->counts[i] = memset (data->counts[i], 0, fastq->size * sizeof (**data->counts));
        }
      data->current_size = fastq->size;
    }
  else if (data->current_size < fastq->size)
    {
      data->ns   = g_realloc (data->ns, fastq->size * sizeof (*data->ns));
      for (i = 0; i < N_LETTERS; i++)
        data->counts[i] = g_realloc (data->counts[i], fastq->size * sizeof (**data->counts));

      for (i = data->current_size; i < fastq->size; i++)
        {
          int j;

          for (j = 0; j < N_LETTERS; j++)
            data->counts[j][i] = 0;

          data->ns[i] = 0;
        }
      data->current_size = fastq->size;
    }

  for (i = 0; i < fastq->size; i++)
    {
      switch (fastq->seq[i])
        {
          case 'A':
          case 'a':
              data->counts[NUC_A][i]++;
              break;
          case 't':
          case 'T':
              data->counts[NUC_T][i]++;
              break;
          case 'g':
          case 'G':
              data->counts[NUC_G][i]++;
              break;
          case 'c':
          case 'C':
              data->counts[NUC_C][i]++;
              break;
          default:
              data->counts[NUC_N][i]++;
              break;
        }
      data->ns[i]++;
    }
  return 1;
}

static void
print_letters (CallbackData *data)
{
  int i;

  g_print ("A\tT\tG\tC\tN\n");
  for (i = 0; i < data->current_size; i++)
    g_print ("%f\t%f\t%f\t%f\t%f\n",
             ((float)data->counts[NUC_A][i]) / data->ns[i],
             ((float)data->counts[NUC_T][i]) / data->ns[i],
             ((float)data->counts[NUC_G][i]) / data->ns[i],
             ((float)data->counts[NUC_C][i]) / data->ns[i],
             ((float)data->counts[NUC_N][i]) / data->ns[i]);

  /* Free stuff up */
  for (i = 0; i < N_LETTERS; i++)
    g_free (data->counts[i]);
  g_free (data->ns);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
