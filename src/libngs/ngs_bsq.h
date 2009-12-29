/**
 *
 */

#ifndef __NGS_BSQ_H__
#define __NGS_BSQ_H__

#include <glib.h>

/************/
/* FastqSeq */
/************/

typedef enum
{
  BSQ_STRAND_W = 0,
  BSQ_STRAND_WC,
  BSQ_STRAND_C,
  BSQ_STRAND_CC,
  BSQ_STRAND_NB,
}
BsqStrand;

typedef enum
{
  BSQ_MAP_UM = 0,
  BSQ_MAP_MA,
  BSQ_MAP_OF,
  BSQ_MAP_NM,
  BSQ_MAP_QC,
  BSQ_MAP_FLAG_NB
}
BsqMapFlag;

extern char *strand_names[];
extern char *map_flag_names[];

typedef struct _BsqRecord BsqRecord;

struct _BsqRecord
{
  char      *name;
  char      *seq;
  char      *ref;
  char      *mis_info;
  char      *mC_loc;
  long int   loc;
  int        size;
  int        n_mis;
  BsqStrand  strand;
  BsqMapFlag flag;
};

BsqRecord*  bsq_record_new    (void);

void        bsq_record_free   (BsqRecord *rec);

/**
 * Signature for the functions to be called by iter_bsq.
 * If the function returns 0, the iteration is interupted,
 * otherwise, the iteration continues.
 */

typedef int (*BsqIterFunc) (BsqRecord *rec,
                            void      *data);

/**
 * Iterates over all FastqSeq in a file.
 * If path is '-', reads from stdin.
 */

void iter_bsq (char         *path,
               BsqIterFunc   func,
               void         *data,
               GError      **error);

/**
 * Get the option group for the fastq parsing system
 */

GOptionGroup* get_bsq_option_group (void);


#endif /* __NGS_BSQ_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
