/**
 *
 */

#ifndef __NGS_SEQ_DB_H__
#define __NGS_SEQ_DB_H__

#include <glib.h>


/****************/
/* SeqDBElement */
/****************/

typedef struct _SeqDBElement SeqDBElement;

struct _SeqDBElement
{
  char    *name;
  guint64  offset;
  guint32  size;
};

SeqDBElement *seq_dbelement_new  (void);
void          seq_dbelement_free (SeqDBElement *elem);

/*********/
/* SeqDB */
/*********/

typedef struct _SeqDB SeqDB;

struct _SeqDB
{
  GHashTable   *index;

  gchar        *seqs;
  gchar        *quals;

  unsigned long alloc_size;
  unsigned long alloc_inc;
  unsigned long total_size;

  unsigned int  n_seqs;
};

SeqDB* seq_db_new        (void);
void   seq_db_free       (SeqDB       *db);
void   seq_db_load_fasta (SeqDB       *db,
                          const char  *path,
                          GError     **error);
void   seq_db_load_fastq (SeqDB       *db,
                          const char  *path,
                          GError     **error);

#endif /* __NGS_SEQ_DB_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
