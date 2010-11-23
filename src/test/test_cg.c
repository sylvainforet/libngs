/* Copyright (C) 2010  Sylvain FORET
 *
 * This file is part of libngs.
 *
 * libngs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *                                                                       
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *                                                                       
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
