/**
 *
 */

#include <stdlib.h>
#include <string.h>

#include "bsq.h"

typedef struct  _BsqSummaryData BsqSummaryData;

struct _BsqSummaryData
{
  unsigned long int  counts[2][BSQ_STRAND_NB][BSQ_MAP_FLAG_NB];
  char              *input_path;
  int                mate_pairs;
  int                strands;
};

static int  iter_func     (BsqRecord         *rec,
                           BsqSummaryData    *data);

static void parse_args    (BsqSummaryData    *data,
                           int               *argc,
                           char            ***argv);

static void print_symmary (BsqSummaryData    *data);

static void init_counts   (BsqSummaryData    *data);

int
main (int    argc,
      char **argv)
{
  BsqSummaryData  data;
  GError         *error = NULL;

  parse_args (&data, &argc, &argv);
  init_counts (&data);

  iter_bsq   (data.input_path,
              (BsqIterFunc)iter_func,
              &data,
              &error);
  if (error)
    {
      g_printerr ("[ERROR] Iterating bsq records failed: %s\n", error->message);
      g_error_free (error);
    }
  else
    print_symmary (&data);

  return 0;
}

static void
parse_args (BsqSummaryData    *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      {"mate_pairs", 'p', 0, G_OPTION_ARG_NONE,   &data->mate_pairs, "Count mate pairs separately", NULL},
      {"strands"   , 's', 0, G_OPTION_ARG_NONE,   &data->strands,    "Counts strands separatedly",  NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->mate_pairs = 0;
  data->strands    = 0;

  context = g_option_context_new ("FILE - Creates a summary of the mapping results");
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
      g_printerr ("[ERROR] No input file provided\n");
      exit (1);
    }
  data->input_path = (*argv)[1];
}

static int
iter_func (BsqRecord      *rec,
           BsqSummaryData *data)
{
#if 0
  if (rec->flag == BSQ_MAP_UM ||
      rec->flag == BSQ_MAP_MA ||
      rec->flag == BSQ_MAP_OF)
    g_print ("%s\t"  /* name */
             "%s\t"  /* seq */
             "%s\t"  /* flag */
             "%s\t"  /* ref */
             "%ld\t" /* loc */
             "%s\t"  /* strand */
             "%d\t"  /* mismatches */
             "%s\t"  /* mismatch_info */
             "%s\n", /* mC_loc */
             rec->name,
             rec->seq,
             map_flag_names[rec->flag],
             rec->ref,
             rec->loc,
             strand_names[rec->strand],
             rec->n_mis,
             rec->mis_info,
             rec->mC_loc);
  else
    g_print ("%s\t"  /* name */
             "%s\t"  /* seq */
             "%s\n", /* flag */
             rec->name,
             rec->seq,
             map_flag_names[rec->flag]);

#else

  int mate_pair_idx = 0;
  int strand_idx    = 0;

  if (data->mate_pairs)
    {
      int name_len = strlen (rec->name);

      if (name_len > 1
          && rec->name[name_len - 1] == '2')
        mate_pair_idx = 1;
    }
  if (data->strands)
    if (rec->flag == BSQ_MAP_UM ||
        rec->flag == BSQ_MAP_MA ||
        rec->flag == BSQ_MAP_OF)
    strand_idx = rec->strand;

  data->counts[mate_pair_idx][strand_idx][rec->flag]++;
#endif

  return 1;
}

static void
print_symmary (BsqSummaryData *data)
{
  int i;
  int j;
  int k;

  if (data->mate_pairs)
    {
      for (i = 0; i < 2; i++)
        {
          g_print ("* Counts for mate pairs %d:\n", i + 1);
          if (data->strands)
            {
              for (j = 0; j < BSQ_STRAND_NB; j++)
                {
                  g_print ("  + Counts for strand %s:\n", strand_names[j]);
                  for (k = 0; k < BSQ_MAP_FLAG_NB; k++)
                    g_print ("    %s: %ld\n",
                             map_flag_names[k],
                             data->counts[i][j][k]);
                }
            }
          else
            {
              g_print ("  + Counts for all strand\n");
              for (k = 0; k < BSQ_MAP_FLAG_NB; k++)
                g_print ("    %s: %ld\n",
                         map_flag_names[k],
                         data->counts[i][0][k]);
            }
        }
    }
  else
    {
      g_print ("* Counts for all mate pairs\n");
      if (data->strands)
        {
          for (j = 0; j < BSQ_STRAND_NB; j++)
            {
              g_print ("  + Counts for strand %s:\n", strand_names[j]);
              for (k = 0; k < BSQ_MAP_FLAG_NB; k++)
                g_print ("    %s: %ld\n",
                         map_flag_names[k],
                         data->counts[0][j][k]);
            }
        }
      else
        {
          g_print ("  + Counts for all strand\n");
          for (k = 0; k < BSQ_MAP_FLAG_NB; k++)
            g_print ("    %s: %ld\n",
                     map_flag_names[k],
                     data->counts[0][0][k]);
        }
    }
}

static void
init_counts (BsqSummaryData *data)
{
  int i;
  int j;
  int k;

  for (i = 0; i < 2; i++)
    for (j = 0; j < BSQ_STRAND_NB; j++)
      for (k = 0; k < BSQ_MAP_FLAG_NB; k++)
        data->counts[i][j][k] = 0;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
