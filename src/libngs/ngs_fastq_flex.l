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

struct _FastqIterFlex
{
  yyscan_t      scanner;
  FILE         *file;
  FlexFastqData data;
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
    /*yyextra->fastq->size = MIN (yyleng, yyextra->fastq->size);*/
    if (yyextra->func)
      {
        ret = yyextra->func (yyextra->fastq, yyextra->data);
        fastq_seq_free (yyextra->fastq);
        yyextra->fastq = NULL;
        if (!ret)
            yyterminate ();
      }
    else
      {
        BEGIN (INITIAL);
        return 0;
      }
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
  data.fastq      = NULL;
  yylex_init_extra (&data, &scanner);
  yyset_in (file, scanner);
  yylex (scanner);
  fclose (file);
  yylex_destroy (scanner);
}

FastqIterFlex* fastq_iter_new_flex (const char   *path,
                                    GError      **error)
{
  FastqIterFlex *iter;
  FILE          *file;

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
      return NULL;
    }

  iter              = g_slice_new (FastqIterFlex);
  iter->data.fastq  = NULL;
  iter->data.func   = NULL;
  iter->data.data   = NULL;
  iter->scanner     = NULL;
  iter->file        = file;
  yylex_init_extra (&iter->data, &iter->scanner);
  yyset_in (file, iter->scanner);

  return iter;
}

FastqSeq* fastq_iter_next_flex (FastqIterFlex *iter)
{
  if (!iter->scanner)
    return NULL;
  if (iter->data.fastq != NULL)
    {
      fastq_seq_free (iter->data.fastq);
      iter->data.fastq = NULL;
    }
  yylex (iter->scanner);
  if (iter->data.fastq == NULL)
    {
      fclose (iter->file);
      iter->file = NULL;
      yylex_destroy (iter->scanner);
      iter->scanner = NULL;
    }
  return iter->data.fastq;
}

void fastq_iter_free_flex (FastqIterFlex *iter)
{
  if (iter->data.fastq != NULL)
    {
      fastq_seq_free (iter->data.fastq);
      iter->data.fastq = NULL;
    }
  if (iter->file != NULL)
    {
      fclose (iter->file);
      iter->file = NULL;
      yylex_destroy (iter->scanner);
    }
  g_slice_free (FastqIterFlex, iter);
}

/* vim:expandtab:ts=4:sw=4:
*/
