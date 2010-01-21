/**
 *
 */

#include <stdlib.h>

#include "ngs_methylation.h"

int
main (int    argc,
      char **argv)
{
  GError        *error = NULL;
  GTimer        *timer;
  SeqDB         *ref;
  RefMethCounts *counts;

  if (argc < 4)
    {
      g_printerr ("Usage: test_cg REF IN OUT\n");
      exit (1);
    }

  timer = g_timer_new ();
  ref = seq_db_new ();
  seq_db_load_fasta (ref, argv[1], &error);
  g_timer_stop (timer);

  if (error)
    {
      g_printerr ("[ERROR] Failed to load reference `%s': %s\n",
                  argv[1],
                  error->message);
      exit (1);
    }
  g_print ("Loaded reference in %.3f sec\n", g_timer_elapsed (timer, NULL));

  g_timer_start (timer);
  counts = ref_meth_counts_load (ref, argv[2], &error);
  g_timer_stop (timer);

  if (error)
    {
      g_printerr ("[ERROR] Failed to load cg file `%s': %s\n",
                  argv[2],
                  error->message);
      exit (1);
    }
  g_print ("Loaded meth file in %.3f sec\n", g_timer_elapsed (timer, NULL));

  g_timer_start (timer);
  ref_meth_counts_write (counts, ref, argv[3], 0, 0, &error);
  g_timer_stop (timer);

  if (error)
    {
      g_printerr ("[ERROR] Failed to write cg file `%s': %s\n",
                  argv[3],
                  error->message);
      exit (1);
    }

  g_print ("Wrote meth file in %.3f sec\n", g_timer_elapsed (timer, NULL));

  g_timer_destroy (timer);

  return 0;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
