%option fast
%option warn
%option yylineno

%option nounput
%option noinput
%option noyywrap

%option reentrant
%option prefix="FlexFastq_"

%option extra-type="FlexFastqData *"

%{
/**
 * FlexFastq.l
 */

#include "ngs_fastq_flex.h"
#include "ngs_utils.h"

typedef struct _FlexFastqData FlexFastqData;

struct _FlexFastqData
{
  FastqSeq      *fastq;
  FastqIterFunc  func;
  void          *data;
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
    yyextra->fastq       = fastq_seq_new ();
    yyextra->fastq->name = strdup (yytext);
    BEGIN (SEQ);
}

^\+.* {
    BEGIN (QUAL);
}

<SEQ>.* {
    yyextra->fastq->seq  = strdup (yytext);
    yyextra->fastq->size = yyleng;
    BEGIN (INITIAL);
}

<QUAL>.* {
    int ret;
    yyextra->fastq->qual = strdup (yytext);
    yyextra->fastq->size = MIN (yyleng, yyextra->fastq->size);
    ret                  = yyextra->func (yyextra->fastq, yyextra->data);
    fastq_seq_free (yyextra->fastq);
    if (!ret)
      yyterminate ();
    BEGIN (INITIAL);
}

<*>[ \t\n] {} /* Ignore spaces */

%%

void iter_fastq_flex (const char   *path,
                      FastqIterFunc func,
                      void         *func_data,
                      GError      **error)
{
  yyscan_t      scanner;
  FILE         *file;
  FlexFastqData data;

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
