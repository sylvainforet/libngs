%option fast
%option warn
%option yylineno

%option nounput
%option noinput
%option noyywrap

%option reentrant
%option prefix="FlexFastqUgly_"

%option extra-type="FlexFastqData *"

%{
/**
 * FlexFastqUgly.l
 */

#include "ngs_fastq_flex.h"
#include "ngs_utils.h"

#define MAX_STRING_SIZE 512

typedef struct _FlexFastqData FlexFastqData;

struct _FlexFastqData
{
  FastqSeq      fastq;
  FastqIterFunc func;
  void         *data;
};

%}

%x NAME
%x SEQ
%x QUAL

%%

^@ {
    BEGIN (NAME);
}

<NAME>.* {
    size_t n = MIN (MAX_STRING_SIZE, yyleng);
    memcpy (yyextra->fastq.name, yytext, n);
    yyextra->fastq.name[n] = '\0';
    BEGIN (SEQ);
}

^\+.* {
    BEGIN (QUAL);
}

<SEQ>.* {
    size_t n = MIN (MAX_STRING_SIZE, yyleng);
    memcpy (yyextra->fastq.seq, yytext, n);
    yyextra->fastq.seq[n] = '\0';
    yyextra->fastq.size = n;
    BEGIN (INITIAL);
}

<QUAL>.* {
    size_t n = MIN (MAX_STRING_SIZE, yyleng);
    memcpy (yyextra->fastq.qual, yytext, n);
    yyextra->fastq.qual[n] = '\0';
    yyextra->fastq.size = MIN (n, yyextra->fastq.size);
    if (!yyextra->func (&yyextra->fastq, yyextra->data))
      yyterminate ();
    BEGIN (INITIAL);
}

<*>[ \t\n] {} /* Ignore spaces */

%%

void iter_fastq_flex_ugly (const char   *path,
                           FastqIterFunc func,
                           void         *func_data,
                           GError      **error)
{
  yyscan_t      scanner;
  FILE         *file;
  FlexFastqData data;
  char          name[MAX_STRING_SIZE + 1];
  char          seq[MAX_STRING_SIZE + 1];
  char          qual[MAX_STRING_SIZE + 1];

  if (path[0] == '-' && path[1] == '\0')
    file = stdin;
  else
    file = fopen (path, "r");
  if (file == NULL)
    {
      g_set_error (error,
                   NGS_ERROR,
                   NGS_IO_ERROR,
                   "Could not open input file `%s'",
                   path);
      return;
    }

  data.fastq.name = name;
  data.fastq.seq  = seq;
  data.fastq.qual = qual;
  data.func       = func;
  data.data       = func_data;
  yylex_init_extra (&data, &scanner);
  yyset_in (file, scanner);
  yylex (scanner);
  fclose (file);
  yylex_destroy (scanner);
}

/* vim:expandtab:ts=4:sw=4:
*/
