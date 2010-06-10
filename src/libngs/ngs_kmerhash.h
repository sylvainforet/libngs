/**
 * Based on Glib's ghash table
 */

#ifndef __NGS_KMERHASH_H__
#define __NGS_KMERHASH_H__

#include <glib.h>

#include "ngs_memalloc.h"

typedef struct _KmerHashNode KmerHashNode;

struct _KmerHashNode
{
  unsigned char *kmer;
  /* If key_hash == 0, node is not in use
   * If key_hash == 1, node is a tombstone
   * If key_hash >= 2, node contains data */
  gulong  key_hash;
  gulong  count;
};

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

  GHashFunc        hash_func;
  GEqualFunc       key_equal_func;
};

typedef struct _KmerHashTableIter KmerHashTableIter;

struct _KmerHashTableIter
{
  KmerHashTable *hash_table;
  glong          position;
};

KmerHashTable* kmer_hash_table_new                 (GHashFunc            hash_func,
                                                    GEqualFunc           key_equal_func,
                                                    gsize                kmer_bytes);

void           kmer_hash_table_destroy             (KmerHashTable       *hash_table);

void           kmer_hash_table_insert              (KmerHashTable       *hash_table,
                                                    unsigned char       *kmer);

KmerHashNode*  kmer_hash_table_lookup              (KmerHashTable       *hash_table,
                                                    const unsigned char *kmer);

void           kmer_hash_table_iter_init           (KmerHashTableIter   *iter,
                                                    KmerHashTable       *hash_table);

KmerHashNode*  kmer_hash_table_iter_next           (KmerHashTableIter   *iter);

#endif /* __NGS_KMERHASH_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
