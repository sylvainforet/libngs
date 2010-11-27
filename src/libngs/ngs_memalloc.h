/**
 *
 */

#ifndef __NGS_MEMALLOC_H__
#define __NGS_MEMALLOC_H__

#include <glib.h>

/**
 * A memory allocator for memmory that does not need to be freed (NF: No Free),
 * except at the very end.
 */

typedef struct _MemAllocNF MemAllocNF;

struct _MemAllocNF
{
  GSList *blocks;

  gsize   elem_size;
  gsize   elem_per_block;
  gsize   block_size;

  gsize   current;
};

MemAllocNF* memallocnf_new  (gsize       elem_size,
                             gsize       elem_per_block);

void        memallocnf_free (MemAllocNF *mem);

void*       memallocnf_add  (MemAllocNF *mem,
                             const void *data);

void*       memallocnf_get  (MemAllocNF *mem);

#endif /* __NGS_MEMALLOC_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
