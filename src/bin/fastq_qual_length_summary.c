/**
 *
 */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ngs_fastq.h"


#define N_QUAL 40

typedef struct _CallbackData CallbackData;

struct _CallbackData
{
  unsigned long int **qual;
  unsigned long int  *ns;

  char               *input_path;

  unsigned int        current_size;
};

static void parse_args (CallbackData      *data,
                        int               *argc,
                        char            ***argv);

static int  iter_func  (FastqSeq          *fastq,
                        CallbackData      *data);

static void print_qual (CallbackData      *data);

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

  print_qual (&data);

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

  context = g_option_context_new ("FILE - Computes distribution of the base qualities depending on the position");
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
  data->qual         = NULL;
  data->ns           = NULL;
}

static int
iter_func (FastqSeq     *fastq,
           CallbackData *data)
{
  int i;

  if (data->current_size == 0)
    {
      data->qual = g_malloc (fastq->size * sizeof (*data->qual));
      data->ns   = g_malloc (fastq->size * sizeof (*data->ns));
      data->ns   = memset (data->ns, 0, fastq->size * sizeof (*data->ns));
      for (i = 0; i < fastq->size; i++)
        {
          data->qual[i] = g_malloc (N_QUAL * sizeof (**data->qual));
          data->qual[i] = memset (data->qual[i], 0, N_QUAL * sizeof (**data->qual));
        }
      data->current_size = fastq->size;
    }
  else if (data->current_size < fastq->size)
    {
      data->qual = g_realloc (data->qual, fastq->size * sizeof (*data->qual));
      data->ns   = g_realloc (data->ns, fastq->size * sizeof (*data->ns));
      for (i = data->current_size; i < fastq->size; i++)
        {
          data->qual[i] = g_malloc (N_QUAL * sizeof (**data->qual));
          data->qual[i] = memset (data->qual[i], 0, N_QUAL * sizeof (**data->qual));
          data->ns[i]   = 0;
        }
      data->current_size = fastq->size;
    }

  for (i = 0; i < fastq->size; i++)
    {
      data->qual[i][fastq->qual[i] - FASTQ_QUAL_0]++;
      data->ns[i]++;
    }

  return 1;
}

static void
print_qual (CallbackData *data)
{
  int i;
  int j;

  for (i = 0; i < data->current_size; i++)
    {
      for (j = 0; j < N_QUAL; j++)
        {
          if (j > 0)
            g_print (" ");
          g_print ("%f", ((double)data->qual[i][j]) / data->ns[i]);
        }
      g_print ("\n");
    }

  /* Free stuff up */
  for (i = 0; i < data->current_size; i++)
    g_free (data->qual[i]);
  g_free (data->qual);
  g_free (data->ns);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
