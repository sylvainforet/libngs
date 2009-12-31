/**
 *
 */

#include <stdlib.h>
#include <string.h>

#include "ngs_bsq.h"


typedef struct _Bsq2GffData Bsq2GffData;

struct _Bsq2GffData
{
  char *source;
  char *input_path;
};

static void parse_args (Bsq2GffData       *data,
                        int               *argc,
                        char            ***argv);

static int  iter_func  (BsqRecord         *rec,
                        Bsq2GffData       *data);

int
main (int    argc,
      char **argv)
{
  Bsq2GffData  data;
  GError      *error = NULL;

  parse_args (&data, &argc, &argv);

  iter_bsq (data.input_path,
            (BsqIterFunc)iter_func,
            &data,
            &error);

  if (error)
    {
      g_printerr ("[ERROR] Iterating sequences failed: %s\n", error->message);
      g_error_free (error);
      error = NULL;
      return 1;
    }

  return 0;
}

static void
parse_args (Bsq2GffData    *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      /* TODO Add an option to select what to print depending on the BsqMapFlag */
      /* TODO Add an option to select what to print depending on the Strand */
      {"source", 's', 0, G_OPTION_ARG_STRING, &data->source, "Read source", NULL},
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  data->source = NULL;

  context = g_option_context_new ("FILE - Converts a bsq file into a gff file");
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

  if (!data->source)
    data->source = "bsread";
}

static int
iter_func (BsqRecord   *rec,
           Bsq2GffData *data)
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
      g_print ("%s\t"     /* Reference */
               "%s\t"     /* Source */
               "read\t"   /* Method */
               "%ld\t"    /* From */
               "%ld\t"    /* To */
               ".\t"      /* Score */
               "%s\t"     /* Strand */
               ".\t"      /* Phase */
               "read %s\n",    /* Group */
               rec->ref,
               data->source,
               rec->loc,
               rec->loc + rec->size,
               (rec->strand == BSQ_STRAND_W || rec->strand == BSQ_STRAND_C) ? "+" : "-",
               rec->name);
    }
  return 1;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
