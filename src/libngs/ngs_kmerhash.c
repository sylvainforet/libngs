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
 * Based on Glib's ghash table
 */

#include <string.h>

#include "ngs_kmerhash.h"

#define HASH_TABLE_MIN_SHIFT 3  /* 1 << 3 == 8 buckets */

/* Each table size has an associated prime modulo (the first prime
 * lower than the table size) used to find the initial bucket. Probing
 * then works modulo 2^n. The prime modulo is necessary to get a
 * good distribution with poor hash functions. */
static const glong prime_mod[] =
{
  1ul,          /* For 1 << 0 */
  2ul,
  3ul,
  7ul,
  13ul,
  31ul,
  61ul,
  127ul,
  251ul,
  509ul,
  1021ul,
  2039ul,
  4093ul,
  8191ul,
  16381ul,
  32749ul,
  65521ul,      /* For 1 << 16 */
  131071ul,
  262139ul,
  524287ul,
  1048573ul,
  2097143ul,
  4194301ul,
  8388593ul,
  16777213ul,
  33554393ul,
  67108859ul,
  134217689ul,
  268435399ul,
  536870909ul,
  1073741789ul,
  2147483647ul,  /* For 1 << 31 */
  4294967291ul,
  8589934583ul,
  17179869143ul,
  34359738337ul,
  68719476731ul,
  137438953447ul,
  274877906899ul,
  549755813881ul,
  1099511627689ul,
  2199023255531ul,
  4398046511093ul,
  8796093022151ul,
  17592186044399ul,
  35184372088777ul,
  70368744177643ul,
  140737488355213ul,
  281474976710597ul,
  562949953421231ul,
  1125899906842597ul,
  2251799813685119ul,
  4503599627370449ul,
  9007199254740881ul,
  18014398509481951ul,
  36028797018963913ul,
  72057594037927931ul,
  144115188075855859ul,
  288230376151711717ul,
  576460752303423433ul,
  1152921504606846883ul,
  2305843009213693951ul,
  4611686018427387847ul,
  9223372036854775783ul,
  18446744073709551557ul
};

static void
kmer_hash_table_set_shift (KmerHashTable *hash_table,
                           gint           shift)
{
  gint i;
  gulong mask = 0;

  hash_table->size = 1 << shift;
  hash_table->mod  = prime_mod[shift];

  for (i = 0; i < shift; i++)
    {
      mask <<= 1;
      mask |= 1;
    }

  hash_table->mask = mask;
}

static gint
kmer_hash_table_find_closest_shift (glong n)
{
  gint i;

  for (i = 0; n; i++)
    n >>= 1;

  return i;
}

static void
kmer_hash_table_set_shift_from_size (KmerHashTable *hash_table,
                                     glong          size)
{
  gint shift;

  shift = kmer_hash_table_find_closest_shift (size);
  shift = MAX (shift, HASH_TABLE_MIN_SHIFT);

  kmer_hash_table_set_shift (hash_table, shift);
}

KmerHashNode*
kmer_hash_table_lookup (KmerHashTable       *hash_table,
                        const unsigned char *kmer)
{
  KmerHashNode *node;
  gulong        node_index;
  gulong        hash_value;
  gulong        step = 0;

  /* Empty buckets have hash_value set to 0, and for tombstones, it's 1.
   * We need to make sure our hash value is not one of these. */
  hash_value = (* hash_table->hash_func) (kmer, hash_table->kmer_bytes);
  if (G_UNLIKELY (hash_value <= 1))
    hash_value = 2;

  node_index = hash_value % hash_table->mod;
  node       = &hash_table->nodes[node_index];

  while (node->key_hash)
    {
      /*  We first check if our full hash values
       *  are equal so we can avoid calling the full-blown
       *  key equality function in most cases.
       */
      if (node->key_hash == hash_value)
        if (hash_table->key_equal_func (kmer, &node->kmer, hash_table->kmer_bytes))
          break;
      step++;
      node_index += step;
      node_index &= hash_table->mask;
      node        = &hash_table->nodes[node_index];
    }
  return node->key_hash ? node : NULL;
}


static inline gulong
kmer_hash_table_lookup_node_for_insertion (KmerHashTable       *hash_table,
                                           const unsigned char *kmer,
                                           gulong              *hash_return)
{
  KmerHashNode *node;
  gulong        node_index;
  gulong        hash_value;
  gulong        first_tombstone = 0;
  gboolean      have_tombstone  = FALSE;
  gulong        step            = 0;

  /* Empty buckets have hash_value set to 0, and for tombstones, it's 1.
   * We need to make sure our hash value is not one of these. */
  hash_value = (* hash_table->hash_func) (kmer, hash_table->kmer_bytes);
  if (G_UNLIKELY (hash_value <= 1))
    hash_value = 2;

  *hash_return = hash_value;
  node_index   = hash_value % hash_table->mod;
  node         = &hash_table->nodes[node_index];

  while (node->key_hash)
    {
      /*  We first check if our full hash values
       *  are equal so we can avoid calling the full-blown
       *  key equality function in most cases.
       */
      if (node->key_hash == hash_value)
        {
          if (hash_table->key_equal_func (kmer, &node->kmer, hash_table->kmer_bytes))
            return node_index;
        }
      else if (node->key_hash == 1 && !have_tombstone)
        {
          first_tombstone = node_index;
          have_tombstone  = TRUE;
        }
      step++;
      node_index += step;
      node_index &= hash_table->mask;
      node        = &hash_table->nodes[node_index];
    }
  if (have_tombstone)
    return first_tombstone;

  return node_index;
}

static void
kmer_hash_table_resize (KmerHashTable *hash_table)
{
/*
  KmerHashNode *new_nodes;
  glong old_size;
  glong i;

  old_size = hash_table->size;
  kmer_hash_table_set_shift_from_size (hash_table, hash_table->nnodes * 2);

  new_nodes = g_new0 (KmerHashNode, hash_table->size);

  for (i = 0; i < old_size; i++)
    {
      KmerHashNode *node = &hash_table->nodes[i];
      KmerHashNode *new_node;
      gulong        hash_val;
      gulong        step = 0;

      if (node->key_hash <= 1)
        continue;

      hash_val = node->key_hash % hash_table->mod;
      new_node = &new_nodes[hash_val];

      while (new_node->key_hash)
        {
          step++;
          hash_val += step;
          hash_val &= hash_table->mask;
          new_node  = &new_nodes[hash_val];
        }

      *new_node = *node;
    }

  g_free (hash_table->nodes);
  hash_table->nodes = new_nodes;
  hash_table->noccupied = hash_table->nnodes;
*/

#define IS_TOUCHED(i) ((touched[(i) / 8] >> ((i) & 7)) & 1)
#define SET_TOUCHED(i) (touched[(i) / 8] |= (1 << ((i) & 7)))

  const gulong   old_size = hash_table->size;
  unsigned char *touched;
  glong          i;

  kmer_hash_table_set_shift_from_size (hash_table, hash_table->nnodes * 2);
  hash_table->nodes = g_realloc (hash_table->nodes, hash_table->size * sizeof (*hash_table->nodes));
  memset (hash_table->nodes + old_size, 0, (hash_table->size - old_size) * sizeof (*hash_table->nodes));
  touched = g_malloc0 ((old_size / 8 + 1) * sizeof (*touched));

  for (i = 0; i < old_size; i++)
    {
      KmerHashNode *node;
      KmerHashNode *next_node;
      KmerHashNode tmp_node;
      gulong       hash_val;
      glong        step;

      if (IS_TOUCHED (i))
        continue;
      node = &hash_table->nodes[i];
      if (node->key_hash == 0)
        continue;
      hash_val = node->key_hash % hash_table->mod;
      if (hash_val == i)
        {
          SET_TOUCHED (i);
          continue;
        }
      tmp_node = *node;
      next_node = &hash_table->nodes[hash_val];
      node->key_hash = 0;
      step = 0;
      while (next_node->key_hash)
        {
          /* Hit an untouched element */
          if (hash_val < old_size && !IS_TOUCHED (hash_val))
            {
              KmerHashNode tmpp;

              tmpp = tmp_node;
              tmp_node = *next_node;
              *next_node = tmpp;
              SET_TOUCHED (hash_val);
              hash_val = tmp_node.key_hash % hash_table->mod;
              step = 0;
            }
          else
            {
              step++;
              hash_val += step;
              hash_val &= hash_table->mask;
            }
          next_node = &hash_table->nodes[hash_val];
        }
      *next_node = tmp_node;
    }
  g_free (touched);
#undef IS_TOUCHED
#undef SET_TOUCHED
}

static inline void
kmer_hash_table_maybe_resize (KmerHashTable *hash_table)
{
  glong noccupied = hash_table->noccupied;
  glong size = hash_table->size;

  if ((size > hash_table->nnodes * 4 && size > 1 << HASH_TABLE_MIN_SHIFT) ||
      (size <= noccupied + (noccupied / 16)))
    kmer_hash_table_resize (hash_table);
}

KmerHashTable*
kmer_hash_table_new (gsize kmer_bytes)
{
  KmerHashTable *hash_table;

  if (kmer_bytes <= KMER_VAL_BYTES)
    hash_table = kmer_hash_table_new_full (kmer_hash_32bp,
                                           kmer_equal_32bp,
                                           kmer_bytes);
  else
    hash_table = kmer_hash_table_new_full (kmer_hash_generic,
                                           kmer_equal_generic,
                                           kmer_bytes);

  return hash_table;
}

KmerHashTable*
kmer_hash_table_new_full (KmerHashFunc  hash_func,
                          KmerEqualFunc key_equal_func,
                          gsize         kmer_bytes)
{
  KmerHashTable *hash_table;

  hash_table                     = g_slice_new (KmerHashTable);
  kmer_hash_table_set_shift (hash_table, HASH_TABLE_MIN_SHIFT);
  hash_table->nnodes             = 0;
  hash_table->noccupied          = 0;
  hash_table->hash_func          = hash_func;
  hash_table->key_equal_func     = key_equal_func;
  hash_table->nodes              = g_new0 (KmerHashNode, hash_table->size);
  hash_table->kmer_bytes         = kmer_bytes;
  if (kmer_bytes > KMER_VAL_BYTES)
    hash_table->allocator        = memallocnf_new (hash_table->kmer_bytes, 1024 * 1024 * 512);
  else
    hash_table->allocator        = NULL;

  return hash_table;
}

void
kmer_hash_table_iter_init (KmerHashTableIter *iter,
                           KmerHashTable     *hash_table)
{
  iter->hash_table = hash_table;
  iter->position   = -1;
}

KmerHashNode*
kmer_hash_table_iter_next (KmerHashTableIter *iter)
{
  KmerHashNode *node      = NULL;
  glong         position  = iter->position;

  do
    {
      position++;
      if (position >= iter->hash_table->size)
        {
          iter->position = position;
          return NULL;
        }
      node = &iter->hash_table->nodes[position];
    }
  while (node->key_hash <= 1);
  iter->position = position;

  return node;
}

void
kmer_hash_table_destroy (KmerHashTable *hash_table)
{
  if (hash_table)
    {
      if (hash_table->allocator)
        memallocnf_free (hash_table->allocator);
      g_slice_free (KmerHashTable, hash_table);
    }
}

void
kmer_hash_table_insert (KmerHashTable *hash_table,
                        unsigned char *kmer)
{
  KmerHashNode *node;
  gulong node_index;
  gulong key_hash;
  gulong old_hash;

  node_index = kmer_hash_table_lookup_node_for_insertion (hash_table, kmer, &key_hash);
  node       = &hash_table->nodes[node_index];
  old_hash   = node->key_hash;

  if (old_hash > 1)
    node->count++;
  else
    {
      if (hash_table->kmer_bytes > KMER_VAL_BYTES)
        node->kmer.kmer_ptr = memallocnf_add (hash_table->allocator, kmer);
      else
        {
          int i;

          for (i = 4; i < hash_table->kmer_bytes; i += 4)
            {
              const int idx1 = i - 4;
              const int idx2 = i - 3;
              const int idx3 = i - 2;
              const int idx4 = i - 1;

              node->kmer.kmer_val[idx1] = kmer[idx1];
              node->kmer.kmer_val[idx2] = kmer[idx2];
              node->kmer.kmer_val[idx3] = kmer[idx3];
              node->kmer.kmer_val[idx4] = kmer[idx4];
            }
          for (i -= 4; i < hash_table->kmer_bytes; i++)
            node->kmer.kmer_val[i] = kmer[i];
        }
      node->count    = 1;
      node->key_hash = key_hash;
      hash_table->nnodes++;
      if (old_hash == 0)
        {
          /* We replaced an empty node, and not a tombstone */
          hash_table->noccupied++;
          kmer_hash_table_maybe_resize (hash_table);
        }
    }
}

gulong
kmer_hash_generic (const unsigned char *kmer,
                   gsize                size)
{
  unsigned long int hash = 0;
  unsigned int      i    = 0;

  do
    {
      hash <<= 8;
      hash  |= kmer[i];
    }
  while (++i < size);

  return ((glong)(hash * 2654435769U));;
}

int
kmer_equal_generic (const unsigned char *kmer1,
                    KmerHashKmer        *kmer2,
                    gsize                size)
{
  unsigned int i;

  for (i = 0; i < size; i++)
    if (kmer1[i] != kmer2->kmer_ptr[i])
      return 0;

  return 1;
}

gulong
kmer_hash_32bp (const unsigned char *kmer,
                gsize                size)
{
  gulong key = *((gulong*)kmer);

  key = (~key) + (key << 21);
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8);
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4);
  key = key ^ (key >> 28);
  key = key + (key << 31);

  return key;
}

/* Completely unrolled */
int
kmer_equal_32bp (const unsigned char *kmer1,
                 KmerHashKmer        *kmer2,
                 gsize                size)
{
  const int b1 = (kmer1[0] == kmer2->kmer_val[0]);
  const int b2 = (kmer1[1] == kmer2->kmer_val[1]);
  const int b3 = (kmer1[2] == kmer2->kmer_val[2]);
  const int b4 = (kmer1[3] == kmer2->kmer_val[3]);
  const int b5 = (kmer1[4] == kmer2->kmer_val[4]);
  const int b6 = (kmer1[5] == kmer2->kmer_val[5]);
  const int b7 = (kmer1[6] == kmer2->kmer_val[6]);
  const int b8 = (kmer1[7] == kmer2->kmer_val[7]);

  const int b12 = b1 && b2;
  const int b34 = b3 && b4;
  const int b56 = b5 && b6;
  const int b78 = b7 && b8;

  const int b1234 = b12 && b34;
  const int b5678 = b56 && b78;

  return b1234 && b5678;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
