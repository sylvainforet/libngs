%option fast
%option warn
%option yylineno

%option nounput
%option noinput
%option noyywrap

%option reentrant
%option prefix="FlexBsqUgly_"

%option extra-type="FlexBsqData *"

%{
/**
 * FlexBsq.l
 */

#include "ngs_bsq_flex.h"
#include "ngs_utils.h"

#define MAX_STRING_SIZE 512

typedef struct _FlexBsqData FlexBsqData;

struct _FlexBsqData
{
  BsqRecord   rec;
  BsqIterFunc func;
  void       *data;
};

%}

%x SEQ
%x FLAG
%x IGNORE
%x REF
%x LOC
%x STRAND
%x NMIS
%x MISINF
%x MCLOC

%%

^[^ \t\n]+ {
    size_t n = MIN (MAX_STRING_SIZE, yyleng);
    memcpy (yyextra->rec.name, yytext, n);
    yyextra->rec.name[n] = '\0';
    BEGIN (SEQ);
}

<SEQ>[^ \t\n]+ {
    size_t n = MIN (MAX_STRING_SIZE, yyleng);
    memcpy (yyextra->rec.seq, yytext, n);
    yyextra->rec.seq[n] = '\0';
    BEGIN (FLAG);
}

<FLAG>[^ \t\n]+ {
    switch (yytext[0])
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
    BEGIN (REF);
}

<REF>[^ \t\n]+ {
    size_t n = MIN (MAX_STRING_SIZE, yyleng);
    memcpy (yyextra->rec.ref, yytext, n);
    yyextra->rec.ref[n] = '\0';
    BEGIN (LOC);
}

<LOC>[^ \t\n]+ {
    yyextra->rec.loc = atol (yytext);
    BEGIN (STRAND);
}

<STRAND>[^ \t\n]+ {
    if (yytext[0] == '+')
      {
        if (yytext[1] == '+')
          yyextra->rec.strand = BSQ_STRAND_W;
        else
          yyextra->rec.strand = BSQ_STRAND_WC;
      }
    else
      {
        if (yytext[1] == '+')
          yyextra->rec.strand = BSQ_STRAND_C;
        else
          yyextra->rec.strand = BSQ_STRAND_CC;
      }
    BEGIN (NMIS);
}

<NMIS>[^ \t\n]+ {
    yyextra->rec.n_mis = atoi (yytext);
    BEGIN (MISINF);
}

<MISINF>[^ \t\n]+ {
    size_t n = MIN (MAX_STRING_SIZE, yyleng);
    memcpy (yyextra->rec.mis_info, yytext, n);
    yyextra->rec.mis_info[n] = '\0';
    BEGIN (MCLOC);
}

<MCLOC>[^ \t\n]+ {
    size_t n = MIN (MAX_STRING_SIZE, yyleng);
    memcpy (yyextra->rec.mC_loc, yytext, n);
    yyextra->rec.mC_loc[n] = '\0';
    BEGIN (IGNORE);
}

<IGNORE>.* {
}

<*>\t {
    /* Ignore tabs */
}

<*>\n {
    if (!yyextra->func (&yyextra->rec, yyextra->data))
      yyterminate ();
    BEGIN (INITIAL);
}

%%

void iter_bsq_flex_ugly (char         *path,
                         BsqIterFunc   func,
                         void         *func_data,
                         GError      **error)
{
  yyscan_t      scanner;
  FILE         *file;
  FlexBsqData   data;
  char          name[MAX_STRING_SIZE + 1];
  char          seq[MAX_STRING_SIZE + 1];
  char          ref[MAX_STRING_SIZE + 1];
  char          mis_info[MAX_STRING_SIZE + 1];
  char          mC_loc[MAX_STRING_SIZE + 1];

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

  data.rec.name     = name;
  data.rec.seq      = seq;
  data.rec.ref      = ref;
  data.rec.mis_info = mis_info;
  data.rec.mC_loc   = mC_loc;
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
