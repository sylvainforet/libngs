/**
 *
 */

#ifndef __NGS_BINSEQ_H__
#define __NGS_BINSEQ_H__

#include <glib.h>

#define BITS_PER_BYTE   8
#define BITS_PER_NUC    2
#define NUCS_PER_BYTE   (BITS_PER_BYTE / BITS_PER_NUC)

extern const unsigned char revcomp_table[];
extern unsigned char       char_to_bin_table[];
extern char                bin_to_char_table[];

/**********/
/* BinSeq */
/**********/

typedef struct _BinSeq BinSeq;

struct _BinSeq
{
  char              *name;
  unsigned char     *seq;
  unsigned long int  size; /* Size in nucleotides */
};


BinSeq*          bin_seq_new                  (const char          *name,
                                               const char          *seq,
                                               unsigned long int    size);

void             bin_seq_free                 (BinSeq              *bin);

unsigned char*   char_to_bin                  (const char          *seq,
                                               unsigned long int    size);

unsigned char*   char_to_bin_prealloc         (unsigned char       *dest,
                                               const char          *src,
                                               unsigned long int    size);

char*            bin_to_char                  (const unsigned char *seq,
                                               unsigned long int    size);

char*            bin_to_char_prealloc         (char                *dest,
                                               const unsigned char *src,
                                               unsigned long int    size);

/* The following functions are for sequences with a number of nucleotides that
 * is a multiple of 4 */

unsigned char*   bin_revcomp_mult4            (const unsigned char *seq,
                                               unsigned long int    size);

unsigned char*   bin_revcomp_prealloc_mult4   (unsigned char       *dest,
                                               const unsigned char *src,
                                               unsigned long int    size);

unsigned char*   bin_revcomp_inplace_mult4    (unsigned char       *seq,
                                               unsigned long int    size);

#endif /* __NGS_BINSEQ_H__ */

/* vim:ft=c:expandtab:sw=4:ts=4:sts=4:cinoptions={.5s^-2n-2(0:
 */
