/**
 * Based on Glib's ghash table
 */

#ifndef __NGS_KMERHASH_H__
#define __NGS_KMERHASH_H__

#include <glib.h>

#include "ngs_memalloc.h"


/***************************/
/* Generic kmer hash table */
/***************************/

#define KMER_VAL_BYTES 8
#define KMER_VAL_NUCS  (KMER_VAL_BYTES * NUCS_PER_BYTE)

typedef union _KmerHashKmer KmerHashKmer;

union _KmerHashKmer
{
  unsigned char *kmer_ptr;
  unsigned char  kmer_val[KMER_VAL_BYTES];
};

typedef struct _KmerHashNode KmerHashNode;

struct _KmerHashNode
{
  KmerHashKmer kmer;

  /* If key_hash == 0, node is not in use
   * If key_hash == 1, node is a tombstone
   * If key_hash >= 2, node contains data */
  gulong  key_hash;
  union
    {
      gulong  count;
      void   *value;
    };
};

typedef gulong (*KmerHashFunc)  (const unsigned char *kmer,
                                 gsize                size);
typedef int    (*KmerEqualFunc) (const unsigned char *kmer1,
                                 KmerHashKmer        *kmer2,
                                 gsize                size);

typedef struct _KmerHashTable  KmerHashTable;

struct _KmerHashTable
{
  KmerHashNode    *nodes;
  MemAllocNF      *allocator;

  glong            size;
  glong            mod;
  gulong           mask;
  glong            nnodes;
  glong            noccupied;  /* nnodes + tombstones */

  gsize            kmer_bytes;

  KmerHashFunc     hash_func;
  KmerEqualFunc    key_equal_func;
};

typedef struct _KmerHashTableIter KmerHashTableIter;

struct _KmerHashTableIter
{
  KmerHashTable *hash_table;
  glong          position;
};

KmerHashTable* kmer_hash_table_new                 (gsize                kmer_bytes);

KmerHashTable* kmer_hash_table_new_full            (KmerHashFunc         hash_func,
                                                    KmerEqualFunc        key_equal_func,
                                                    gsize                kmer_bytes);

void           kmer_hash_table_destroy             (KmerHashTable       *hash_table);

void           kmer_hash_table_increment           (KmerHashTable       *hash_table,
                                                    const unsigned char *kmer);

KmerHashNode*  kmer_hash_table_lookup_or_create    (KmerHashTable       *hash_table,
                                                    const unsigned char *kmer);

KmerHashNode*  kmer_hash_table_lookup              (KmerHashTable       *hash_table,
                                                    const unsigned char *kmer);

void           kmer_hash_table_iter_init           (KmerHashTableIter   *iter,
                                                    KmerHashTable       *hash_table);

KmerHashNode*  kmer_hash_table_iter_next           (KmerHashTableIter   *iter);

gulong         kmer_hash_generic                   (const unsigned char *kmer,
                                                    gsize                size);

int            kmer_equal_generic                  (const unsigned char *kmer1,
                                                    KmerHashKmer        *kmer2,
                                                    gsize                size);

gulong         kmer_hash_32bp                      (const unsigned char *kmer,
                                                    gsize                size);

int            kmer_equal_32bp                     (const unsigned char *kmer1,
                                                    KmerHashKmer        *kmer2,
                                                    gsize                size);

#endif /* __NGS_KMERHASH_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
