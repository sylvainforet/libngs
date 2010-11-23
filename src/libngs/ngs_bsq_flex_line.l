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
%option prefix="FlexBsqLine_"

%option extra-type="FlexBsqData *"

%{
/**
 * FlexBsqLine.l
 */

#include "ngs_bsq_flex.h"
#include "ngs_utils.h"

typedef struct _FlexBsqData FlexBsqData;

struct _FlexBsqData
{
  BsqRecord    rec;
  BsqIterFunc  func;
  void        *data;
};

%}

%%

^.+\n {
    char *tmp;
    char *tmpp;

    /* Clear the structure */
    yyextra->rec.name     = NULL;
    yyextra->rec.seq      = NULL;
    yyextra->rec.ref      = NULL;
    yyextra->rec.mis_info = NULL;
    yyextra->rec.mC_loc   = NULL;
    yyextra->rec.size     = 0;
    yyextra->rec.loc      = 0;
    yyextra->rec.n_mis    = 0;

    tmp = yytext;
    /* Name */
    yyextra->rec.name = tmp;
    while (*tmp != '\t')
      ++tmp;
    *tmp = '\0';

    /* Seq */
    ++tmp;
    yyextra->rec.seq = tmp;
    while (*tmp != '\t')
      ++tmp;
    *tmp = '\0';
    yyextra->rec.size = strlen (yyextra->rec.seq);

    /* Flag */
    ++tmp;
    switch (*tmp)
      {
        case 'U':
            yyextra->rec.flag = BSQ_MAP_UM;
            break;
        case 'M':
            yyextra->rec.flag = BSQ_MAP_MA;
            break;
        case 'O':
            yyextra->rec.flag = BSQ_MAP_OF;
            break;
        case 'N':
            yyextra->rec.flag = BSQ_MAP_NM;
            break;
        case 'Q':
            yyextra->rec.flag = BSQ_MAP_QC;
            break;
        default:
            break;
      }

    if (yyextra->rec.flag == BSQ_MAP_UM ||
        yyextra->rec.flag == BSQ_MAP_MA ||
        yyextra->rec.flag == BSQ_MAP_OF)
      {
        ++tmp; /* 2nd letter */
        ++tmp; /* \t */

        /* Ref */
        ++tmp;
        yyextra->rec.ref = tmp;
        while (*tmp != '\t')
          ++tmp;
        *tmp = '\0';

        /* Loc */
        ++tmp;
        tmpp = tmp;
        while (*tmp != '\t')
          ++tmp;
        *tmp = '\0';
        yyextra->rec.loc = atol (tmpp);

        /* Strand */
        ++tmp;
        tmpp = tmp;
        ++tmp; /* 2nd letter */
        ++tmp; /* \t */
        *tmp = '\0';
        if (tmpp[0] == '+')
          {
            if (tmpp[1] == '+')
              yyextra->rec.strand = BSQ_STRAND_W;
            else
              yyextra->rec.strand = BSQ_STRAND_WC;
          }
        else
          {
            if (tmpp[1] == '+')
              yyextra->rec.strand = BSQ_STRAND_C;
            else
              yyextra->rec.strand = BSQ_STRAND_CC;
          }

        /* n_mis */
        ++tmp;
        tmpp = tmp;
        while (*tmp != '\t')
          ++tmp;
        *tmp = '\0';
        yyextra->rec.n_mis = atoi (tmpp);

        /* mis_info */
        ++tmp;
        yyextra->rec.mis_info = tmp;
        while (*tmp != '\t')
          ++tmp;
        *tmp = '\0';

        /* mC_loc */
        ++tmp;
        yyextra->rec.mC_loc = tmp;
        while (*tmp != '\t' && *tmp != '\n')
          ++tmp;
        *tmp = '\0';
      }

    if (!yyextra->func (&yyextra->rec, yyextra->data))
      yyterminate ();
}

\n {
    /* Ignore empty lines */
}

%%

void iter_bsq_flex_line (char         *path,
                         BsqIterFunc   func,
                         void         *func_data,
                         GError      **error)
{
  yyscan_t      scanner;
  FILE         *file;
  FlexBsqData   data;

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

  data.func         = func;
  data.data         = func_data;
  yylex_init_extra (&data, &scanner);
  yyset_in (file, scanner);
  yylex (scanner);
  fclose (file);
  yylex_destroy (scanner);
}

/* vim:expandtab:ts=4:sw=4:
*/
