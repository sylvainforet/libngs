%option fast
%option warn
%option yylineno

%option nounput
%option noinput
%option noyywrap

%option reentrant
%option prefix="FlexFasta_"

%option extra-type="FlexFastaData *"

%{
/**
 * FlexFasta.l
 */

#include "ngs_fasta_flex.h"

typedef struct _FlexFastaData FlexFastaData;

struct _FlexFastaData
{
  FastaSeq      *fasta;
  FastaIterFunc  func;
  GString       *buffer;
  void          *data;
};

%}

%x NAME
%x HEADER_IGNORE
%x SEQ

%%

<*>>[ \t]* {
    if (yyextra->fasta)
      {
        int ret;
        yyextra->fasta->seq  = yyextra->buffer->str;
        yyextra->fasta->size = yyextra->buffer->len;
        ret = yyextra->func (yyextra->fasta, yyextra->data);
        g_string_free (yyextra->buffer, TRUE);
        yyextra->fasta->seq  = NULL;
        fasta_seq_free (yyextra->fasta);
        yyextra->fasta       = NULL;
        if (!ret)
          yyterminate ();
      }
    yyextra->fasta  = fasta_seq_new ();
    yyextra->buffer = g_string_sized_new (1024);
    BEGIN (NAME);
}

<NAME>[^ \t\n]+ {
    yyextra->fasta->name = strdup (yytext);
    BEGIN (HEADER_IGNORE);
}

<HEADER_IGNORE>.*\n {
    BEGIN (SEQ);
}

<SEQ>[^ \t\n>]+ {
    yyextra->buffer = g_string_append (yyextra->buffer, yytext);
}

<*>[ \t\n] {
    /* Ignore spaces */
}

<SEQ><<EOF>> {
    if (yyextra->fasta)
      {
        yyextra->fasta->seq  = yyextra->buffer->str;
        yyextra->fasta->size = yyextra->buffer->len;
        yyextra->func (yyextra->fasta, yyextra->data);
        g_string_free (yyextra->buffer, TRUE);
        yyextra->fasta->seq  = NULL;
        fasta_seq_free (yyextra->fasta);
        yyextra->fasta = NULL;
      }
    yyterminate ();
}

%%

void iter_fasta_flex (char         *path,
                      FastaIterFunc func,
                      void         *func_data,
                      GError      **error)
{
  yyscan_t      scanner;
  FILE         *file;
  FlexFastaData data;

  if (path[0] == '-' && path[1] == '\0')
    file = stdin;
  else
    {
      file = fopen (path, "r");
      if (file == NULL)
        {
          g_printerr ("[ERROR] Could not open %s\n", path);
          /* TODO Should raise an error instead => have to create my own error domain ... */
          exit (1);
        }
    }

  data.fasta      = NULL;
  data.buffer     = NULL;
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
