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

#include <string.h>

#include "ngs_fasta.h"
#include "ngs_fastq.h"
#include "ngs_seq_db.h"
#include "ngs_utils.h"


static int  iter_load_db_fastq (FastqSeq *fastq,
                                SeqDB    *db);

static int  iter_load_db_fasta (FastaSeq *fasta,
                                SeqDB    *db);

SeqDBElement*
seq_db_element_new (void)
{
  SeqDBElement *elem;

  elem = g_slice_new0 (SeqDBElement);

  return elem;
}

void
seq_db_element_free (SeqDBElement *elem)
{
  if (elem)
    {
      if (elem->name)
        g_free (elem->name);
      g_slice_free (SeqDBElement, elem);
    }
}

#define DEFAULT_ALLOC_INC (16 * 1024 * 1024) /* 2^4 * 2^10 * 2^10 = 2^22 */
SeqDB*
seq_db_new (void)
{
  SeqDB *db;

  db             = g_slice_new (SeqDB);
  db->index      = g_hash_table_new_full (g_str_hash,
                                          g_str_equal,
                                          NULL,
                                          (GDestroyNotify)seq_db_element_free);
  db->seqs       = NULL;
  db->quals      = NULL;
  db->alloc_size = 0;
  db->total_size = 0;
  db->alloc_inc  = DEFAULT_ALLOC_INC;
  db->n_seqs     = 0;

  return db;
}
#undef DEFAULT_ALLOC_INC

void
seq_db_free (SeqDB *db)
{
  if (db)
    {
      if (db->index)
        g_hash_table_destroy (db->index);
      if (db->seqs)
        g_free (db->seqs);
      if (db->quals)
        g_free (db->quals);
      g_slice_free (SeqDB, db);
    }
}

void
seq_db_load_fasta (SeqDB       *db,
                   const char  *path,
                   GError     **error)
{
  iter_fasta (path,
              (FastaIterFunc)iter_load_db_fasta,
              db,
              error);
}

void
seq_db_load_fastq (SeqDB       *db,
                   const char  *path,
                   GError     **error)
{
  iter_fastq (path,
              (FastqIterFunc)iter_load_db_fastq,
              db,
              error);
}

static int
iter_load_db_fastq (FastqSeq *fastq,
                    SeqDB    *db)
{
  SeqDBElement *elem;

  db->n_seqs++;
  if (db->total_size >= db->alloc_size)
    {
      db->alloc_size  = ((db->alloc_size + fastq->size + db->alloc_inc - 1) / db->alloc_inc) * db->alloc_inc;
      db->seqs        = g_realloc (db->seqs,
                                   db->alloc_size * sizeof (*db->seqs));
      db->quals       = g_realloc (db->quals,
                                   db->alloc_size * sizeof (*db->quals));
    }
  elem            = seq_db_element_new ();
  elem->name      = strdup (fastq->name);
  elem->offset    = db->total_size;
  elem->size      = fastq->size;
  memcpy (db->seqs + db->total_size,
          fastq->seq,
          fastq->size * sizeof (*fastq->seq));
  memcpy (db->quals + db->total_size,
          fastq->qual,
          fastq->size * sizeof (*fastq->seq));
  db->total_size += elem->size;
  g_hash_table_insert (db->index,
                       elem->name,
                       elem);

  return 1;
}

static int
iter_load_db_fasta (FastaSeq *fasta,
                    SeqDB    *db)
{
  SeqDBElement *elem;

  db->n_seqs++;
  if (db->total_size + fasta->size >= db->alloc_size)
    {
      db->alloc_size  = ((db->alloc_size + fasta->size + db->alloc_inc - 1) / db->alloc_inc) * db->alloc_inc;
      db->seqs        = g_realloc (db->seqs,
                                   db->alloc_size * sizeof (*db->seqs));
    }
  elem            = seq_db_element_new ();
  elem->name      = strdup (fasta->name);
  elem->offset    = db->total_size;
  elem->size      = fasta->size;
  memcpy (db->seqs + db->total_size,
          fasta->seq,
          fasta->size * sizeof (*fasta->seq));
  db->total_size += elem->size;
  g_hash_table_insert (db->index,
                       elem->name,
                       elem);

  return 1;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
