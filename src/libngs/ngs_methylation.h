/**
 *
 */

#ifndef __NGS_METHYLATION_H__
#define __NGS_METHYLATION_H__

#include "ngs_seq_db.h"


/*************/
/* MethCount */
/*************/

typedef struct _MethCount MethCount;

struct _MethCount
{
  unsigned int n_meth;
  unsigned int n_unmeth;
};


typedef struct _RefMethCounts RefMethCounts;

struct _RefMethCounts
{
  MethCount  *meth_data;
  MethCount **meth_index;
};

RefMethCounts* ref_meth_counts_create        (SeqDB         *ref);

RefMethCounts* ref_meth_counts_load          (SeqDB         *ref,
                                              const char    *path,
                                              GError       **error);

void           ref_meth_counts_add_path      (RefMethCounts *counts,
                                              SeqDB         *ref,
                                              const char    *path,
                                              GError       **error);

void           ref_meth_counts_write         (RefMethCounts *counts,
                                              SeqDB         *ref,
                                              const char    *path,
                                              int            print_letter,
                                              int            print_all,
                                              GError       **error);

void           ref_meth_counts_write_segment (RefMethCounts *counts,
                                              SeqDB         *ref,
                                              GIOChannel    *channel,
                                              const char    *name,
                                              unsigned long  from,
                                              unsigned long  to,
                                              int            print_letter,
                                              int            print_all,
                                              GError       **error);

void           ref_meth_counts_destroy  (RefMethCounts *counts);

#endif /* __NGS_METHYLATION_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
