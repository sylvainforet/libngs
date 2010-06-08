/**
 *
 */

#include <string.h>

#include "ngs_memalloc.h"


static gsize nearest_power_of_2 (gsize n);


MemAllocNF*
memallocnf_new (gsize elem_size,
                gsize elem_per_block)
{
  MemAllocNF *mem;

  mem                 = g_slice_new (MemAllocNF);
  /* TODO Should elem_size be aligned on the 32/64 bits boundary? */
  mem->elem_size      = elem_size;
  mem->block_size     = nearest_power_of_2 (mem->elem_size * elem_per_block);
  mem->elem_per_block = mem->block_size / mem->elem_size;
  mem->current        = 0;
  mem->blocks         = NULL;

  return mem;
}

void
memallocnf_free (MemAllocNF *mem)
{
  if (mem)
    {
      GSList *list;

      for (list = mem->blocks; list; list = list->next)
        g_free (list->data);
      g_slist_free (mem->blocks);
      g_slice_free (MemAllocNF, mem);
    }
}

void*
memallocnf_add (MemAllocNF *mem,
                void       *data)
{
  void *ret;

  if (!mem->blocks || mem->current + mem->elem_size > mem->block_size)
    {
      mem->blocks  = g_slist_prepend (mem->blocks, g_malloc (mem->block_size));
      mem->current = 0;
    }
  ret = memcpy (mem->blocks->data + mem->current, data, mem->elem_size);
  mem->current += mem->elem_size;

  return ret;
}

void*
memallocnf_get (MemAllocNF *mem)
{
  void *ret;

  if (!mem->blocks || mem->current + mem->elem_size > mem->block_size)
    {
      mem->blocks  = g_slist_prepend (mem->blocks, g_malloc (mem->block_size));
      mem->current = 0;
    }

  ret = mem->blocks->data + mem->current;
  mem->current += mem->elem_size;

  return ret;
}

static gsize
nearest_power_of_2 (gsize n)
{
  gsize i = 1;

  while (i < n)
    i <<= 1;

  return i;
}

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
