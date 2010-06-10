/**
 * Based on Glib's ghash table
 */

#include "ngs_kmerhash.h"

#define HASH_TABLE_MIN_SHIFT 3  /* 1 << 3 == 8 buckets */

/* Each table size has an associated prime modulo (the first prime
 * lower than the table size) used to find the initial bucket. Probing
 * then works modulo 2^n. The prime modulo is necessary to get a
 * good distribution with poor hash functions. */
static const glong prime_mod [] =
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
  1099511627689ul
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
kmer_hash_table_find_closest_shift (gint n)
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
  hash_value = (* hash_table->hash_func) (kmer);
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
        if (hash_table->key_equal_func (node->kmer, kmer))
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
  hash_value = (* hash_table->hash_func) (kmer);
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
          if (hash_table->key_equal_func (node->kmer, kmer))
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
  KmerHashNode *new_nodes;
  glong old_size;
  glong i;

  old_size = hash_table->size;
  kmer_hash_table_set_shift_from_size (hash_table, hash_table->nnodes * 2);

  g_printerr (">>>>>> Resizing hash table from % 12ld to % 12ld -- mod: % 12ld, nodes: % 12ld",
          old_size, hash_table->size, hash_table->mod, hash_table->nnodes);
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
  g_printerr (">>> ... DONE\n");
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
kmer_hash_table_new (GHashFunc  hash_func,
                     GEqualFunc key_equal_func,
                     gsize      kmer_bytes)
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
  hash_table->allocator          = memallocnf_new (hash_table->kmer_bytes, 1024 * 1024 * 512);

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
    g_slice_free (KmerHashTable, hash_table);
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
      node->kmer     = memallocnf_add (hash_table->allocator, kmer);
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

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
