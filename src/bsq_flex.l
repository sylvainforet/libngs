%option fast
%option warn
%option yylineno

%option nounput
%option noinput
%option noyywrap

%option reentrant
%option prefix="FlexBsq_"

%option extra-type="FlexBsqData *"

%{
/**
 * FlexBsq.l
 */

#include "bsq_flex.h"

typedef struct _FlexBsqData FlexBsqData;

struct _FlexBsqData
{
  BsqRecord   *rec;
  BsqIterFunc  func;
  void        *data;
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
    yyextra->rec       = bsq_record_new ();
    yyextra->rec->name = strdup (yytext);
    BEGIN (SEQ);
}

<SEQ>[^ \t\n]+ {
    yyextra->rec->seq = strdup (yytext);
    BEGIN (FLAG);
}

<FLAG>[^ \t\n]+ {
    switch (yytext[0])
      {
        case 'U':
            yyextra->rec->flag = BSQ_MAP_UM;
            break;
        case 'M':
            yyextra->rec->flag = BSQ_MAP_MA;
            break;
        case 'O':
            yyextra->rec->flag = BSQ_MAP_OF;
            break;
        case 'N':
            yyextra->rec->flag = BSQ_MAP_NM;
            break;
        case 'Q':
            yyextra->rec->flag = BSQ_MAP_QC;
            break;
        default:
            break;
      }
    BEGIN (REF);
}

<REF>[^ \t\n]+ {
    yyextra->rec->ref = strdup (yytext);
    BEGIN (LOC);
}

<LOC>[^ \t\n]+ {
    yyextra->rec->loc = atol (yytext);
    BEGIN (STRAND);
}

<STRAND>[^ \t\n]+ {
    if (yytext[0] == '+')
      {
        if (yytext[1] == '+')
          yyextra->rec->strand = BSQ_STRAND_W;
        else
          yyextra->rec->strand = BSQ_STRAND_WC;
      }
    else
      {
        if (yytext[1] == '+')
          yyextra->rec->strand = BSQ_STRAND_C;
        else
          yyextra->rec->strand = BSQ_STRAND_CC;
      }
    BEGIN (NMIS);
}

<NMIS>[^ \t\n]+ {
    yyextra->rec->n_mis = atoi (yytext);
    BEGIN (MISINF);
}

<MISINF>[^ \t\n]+ {
    yyextra->rec->mis_info = strdup (yytext);
    BEGIN (MCLOC);
}

<MCLOC>[^ \t\n]+ {
    yyextra->rec->mC_loc = strdup (yytext);
    BEGIN (IGNORE);
}

<IGNORE>.* {
}

<*>\t {
    /* Ignore tabs */
}

<*>\n {
    int ret;

    ret = yyextra->func (yyextra->rec, yyextra->data);
    if (yyextra->rec)
      bsq_record_free (yyextra->rec);
    yyextra->rec = NULL;
    if (!ret)
      yyterminate ();
    BEGIN (INITIAL);
}

%%

void iter_bsq_flex (char         *path,
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
    {
      file = fopen (path, "r");
      if (file == NULL)
        {
          g_printerr ("[ERROR] Could not open %s\n", path);
          /* TODO Should raise an error instead => have to create my own error domain ... */
          exit (1);
        }
    }

  data.rec          = NULL;
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