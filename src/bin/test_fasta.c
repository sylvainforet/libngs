/**
 *
 */

#include <stdlib.h>

#include "fasta.h"

typedef struct  _TestFastaData TestFastaData;

struct _TestFastaData
{
  char *input_path;
  int   width;
};

static void parse_args (TestFastaData    *data,
                        int               *argc,
                        char            ***argv);

static void iter_fasta_dict_func (char *key,
                                  char *value,
                                  void *data);

int
main (int    argc,
      char **argv)
{
  TestFastaData  data;
  GHashTable    *dict;

  parse_args (&data, &argc, &argv);
  dict = load_fasta_dict (data.input_path);
  g_hash_table_foreach (dict, (GHFunc)iter_fasta_dict_func, NULL);

  return 0;
}

static void
parse_args (TestFastaData    *data,
            int               *argc,
            char            ***argv)
{
  GOptionEntry entries[] =
    {
      /* {"width", 'w', 0, G_OPTION_ARG_INT, &data->width, "Line width", NULL}, */
      {NULL}
    };
  GError         *error = NULL;
  GOptionContext *context;

  context = g_option_context_new ("FILE - Reads a fasta file and spits it back out");
  g_option_context_add_group (context, get_fasta_option_group ());
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

static void
iter_fasta_dict_func (char *key,
                      char *value,
                      void *data)
{
  g_print (">%s\n%s\n",
           key,
           value);
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
