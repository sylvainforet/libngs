/**
 *
 */

#include <stdlib.h>

#include "bsq.h"

typedef struct  _BsqSummaryData BsqSummaryData;

struct _BsqSummaryData
{
  unsigned long int ***counts;
  char                *input_path;
  int                  mate_pairs;
  int                  strands;
};

static int  iter_func  (BsqRecord         *rec,
                        BsqSummaryData    *data);

static void parse_args (BsqSummaryData    *data,
                        int               *argc,
                        char            ***argv);

int
main (int    argc,
      char **argv)
{
  BsqSummaryData  data;
  GError         *error = NULL;

  parse_args (&data, &argc, &argv);

  iter_bsq   (data.input_path,
              (BsqIterFunc)iter_func,
              &data,
              &error);
  if (error)
    {
      g_printerr ("[ERROR] Iterating bsq records failed: %s\n", error->message);
      g_error_free (error);
    }

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

  return 1;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
