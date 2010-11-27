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

#include <glib.h>

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

struct _FastaIterFlex
{
  yyscan_t      scanner;
  FILE         *file;
  FlexFastaData data;
};

%}

%x NAME
%x HEADER_IGNORE
%x SEQ

%%

<*>>[ \t]* {
    if (yyextra->fasta)
      {
        yyextra->fasta->seq  = yyextra->buffer->str;
        yyextra->fasta->size = yyextra->buffer->len;
        if (yyextra->func)
          {
            int ret;
            ret                  = yyextra->func (yyextra->fasta, yyextra->data);
            g_string_truncate (yyextra->buffer, 0);
            yyextra->fasta->seq  = NULL;
            fasta_seq_free (yyextra->fasta);
            yyextra->fasta       = NULL;
            if (!ret)
              yyterminate ();
          }
        else
          {
            BEGIN (NAME);
            return 0;
          }
      }
    BEGIN (NAME);
}

<NAME>[^ \t\n]+ {
    if (!yyextra->func && yyextra->fasta)
      {
        g_string_truncate (yyextra->buffer, 0);
        yyextra->fasta->seq  = NULL;
        fasta_seq_free (yyextra->fasta);
        yyextra->fasta       = NULL;
      }
    yyextra->fasta  = fasta_seq_new ();
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
        if (yyextra->func)
          {
            yyextra->func (yyextra->fasta, yyextra->data);
            g_string_free (yyextra->buffer, TRUE);
            yyextra->fasta->seq  = NULL;
            fasta_seq_free (yyextra->fasta);
            yyextra->fasta = NULL;
          }
        else
          return 0;
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
  data.buffer     = g_string_sized_new (4096);
  data.func       = func;
  data.data       = func_data;
  yylex_init_extra (&data, &scanner);
  yyset_in (file, scanner);
  yylex (scanner);
  fclose (file);
  yylex_destroy (scanner);
}

FastaIterFlex* fasta_iter_new_flex (const char   *path,
                                    GError      **error)
{
  FastaIterFlex *iterator;
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

  iterator              = g_slice_new (FastaIterFlex);
  iterator->data.fasta  = NULL;
  iterator->data.buffer = g_string_sized_new (4096);
  iterator->data.func   = NULL;
  iterator->data.data   = NULL;
  iterator->scanner     = NULL;
  iterator->file        = file;
  yylex_init_extra (&iterator->data, &iterator->scanner);
  yyset_in (file, iterator->scanner);

  return iterator;
}

FastaSeq* fasta_iter_next_flex (FastaIterFlex *iter)
{
  if (iter->data.fasta != NULL)
    {
      g_string_truncate (iter->data.buffer, 0);
      iter->data.fasta->seq  = NULL;
      fasta_seq_free (iter->data.fasta);
      iter->data.fasta       = NULL;
    }
  yylex (iter->scanner);
  if (iter->data.fasta == NULL)
    {
      fclose (iter->file);
      iter->file = NULL;
      yylex_destroy (iter->scanner);
      return NULL;
    }
  return iter->data.fasta;
}

void fasta_iter_free_flex (FastaIterFlex *iter)
{
  if (iter->data.fasta != NULL)
    {
      g_string_free (iter->data.buffer, TRUE);
      iter->data.fasta->seq  = NULL;
      fasta_seq_free (iter->data.fasta);
      iter->data.fasta       = NULL;
    }
  if (iter->file != NULL)
    {
      fclose (iter->file);
      iter->file = NULL;
      yylex_destroy (iter->scanner);
    }
  g_slice_free (FastaIterFlex, iter);
}

/* vim:expandtab:ts=4:sw=4:
*/
