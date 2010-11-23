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
#include "ngs_utils.h"

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

void iter_fasta_flex (const char   *path,
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
