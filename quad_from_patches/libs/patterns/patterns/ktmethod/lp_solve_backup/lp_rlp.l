WS1  [ \t]
WS2  [ \t\r\n]
LT   [A-Za-z]
KR   [A-Za-z0-9_\[\]\{\}/.&#$%~'@^]
DI   [0-9]
INNM {DI}+
NM   (({DI}+\.?{DI}*)|(\.{DI}+))([Ee][-+]?{DI}+)?
VR   {LT}{KR}*(<{KR}+>)?
OPEQ =
OPLE [<]=?
OPGE [>]=?
EOW  $|[*+-;<=>]|{WS2}+
INF  {WS2}*[-+][Ii][Nn][Ff]([Ii][Nn][Ii][Tt]([Ee]|[Yy]))?{EOW}
S_OP [-+ \t\n\r]+
MIN  [mM][iI][nN]([iI][mM][iI][zZsS][eE])?{WS1}*:
MAX  [mM][aA][xX]([iI][mM][iI][zZsS][eE])?{WS1}*:
INT  ^{WS1}*[Ii][Nn][Tt]([Ee][Gg][Ee][Rr])?
BIN  ^{WS1}*[Bb][Ii][Nn]([Aa][Rr][Yy])?
SEC  ^{WS1}*([Ss][Ee][Cc])|([Ss][Ee][Mm][Ii]-[Cc][Oo][Nn][Tt][Ii][Nn][Uu][Oo][Uu][Ss])|([Ss][Ee][Mm][Ii])|([Ss][Ee][Mm][Ii][Ss])
SOS  ^{WS1}*([Ss][Oo][Ss][12]*)|([Ss][Ee][Tt][Ss])
FREE ^{WS1}*[Ff][Rr][Ee][Ee]
LBL  {VR}{WS1}*:
SOSD {LBL}:

/*
   made reentrant with help of
   http://www.usualcoding.eu/post/2007/09/03/Building-a-reentrant-parser-in-C-with-Flex/Bison
*/

/*
   Note that a minimum version of flex is needed to be able to compile this.
   Older version don't know the reentrant code.
   Version 2.5.4 is not enough. Probably at least v2.5.31 is needed. Tested with v2.5.35
*/

/*
** We want the scanner to be reentrant, therefore generate no global variables.
** That what the 'reentrant' option is for.
** 'bison-bridge' is used to create a bison compatible scanner and share yylval
*/

%option reentrant bison-bridge
%option noyywrap
%option yylineno

%start COMMENT
%start LINECOMMENT

%%
<INITIAL>"/*" {
  BEGIN COMMENT;
} /* begin skip comment */

<COMMENT>"*/" {
  BEGIN INITIAL;
} /* end skip comment */

<COMMENT>. {
}

<COMMENT>\n {
}

<COMMENT>\r {
}

<INITIAL>"//" {
  BEGIN LINECOMMENT;
} /* begin skip LINECOMMENT */

<LINECOMMENT>\n {
  BEGIN INITIAL;
} /* end skip LINECOMMENT */

<LINECOMMENT>\r {
  BEGIN INITIAL;
} /* end skip LINECOMMENT */

<LINECOMMENT>. {
}

<INITIAL>{WS2}+ {
}

<INITIAL>"," {
  parse_parm *pp = PARM;

  pp->lineno = yylineno;
  return(COMMA);
}

<INITIAL>{MIN} {
  parse_parm *pp = PARM;

  pp->lineno = yylineno;
  return(MINIMISE);
}

<INITIAL>{MAX} {
  parse_parm *pp = PARM;

  pp->lineno = yylineno;
  return(MAXIMISE);
}

<INITIAL>{INNM} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  pv->f = atof((char *)yytext);
  return(INTCONS);
} /* f contains the last float */

<INITIAL>{NM} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  pv->f = atof((char *)yytext);
  return(CONS);
} /* f contains the last float */

<INITIAL>{INF} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;
  char *ptr, c;

  pp->lineno = yylineno;
  pv->f = DEF_INFINITE;
  pv->Sign = 0;
  ptr = (char *)yytext;
  while (isspace(*ptr)) ptr++;
  if(*ptr == '-')
    pv->Sign = 1;
  if(yyleng > 0) {
    c = yytext[yyleng - 1];
    if(!isalnum(c))
      unput(c);
  }
  return(INF);
} /* f contains the last float */

<INITIAL>{S_OP} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;
  int x;

  pp->lineno = yylineno;
  pv->Sign = 0;
  for(x = 0; x < yyleng; x++)
    if(yytext[x] == '-' || yytext[x] == '+')
      pv->Sign = (pv->Sign == (yytext[x] == '+'));
  return (TOK_SIGN);
  /* Sign is TRUE if the sign-string
     represents a '-'. Otherwise Sign
     is FALSE */
}

<INITIAL>{INT} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  if((!pv->Within_int_decl) && (!pv->Within_sec_decl) && (!pv->Within_sos_decl) && (!pv->Within_free_decl)) {
    pv->Within_int_decl = 1;
    pv->Within_sos_decl1 = FALSE;
  }
  return(SEC_INT);
}

<INITIAL>{BIN} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  if((!pv->Within_int_decl) && (!pv->Within_sec_decl) && (!pv->Within_sos_decl) && (!pv->Within_free_decl)) {
    pv->Within_int_decl = 2;
    pv->Within_sos_decl1 = FALSE;
  }
  return(SEC_BIN);
}

<INITIAL>{SEC} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  if((!pv->Within_int_decl) && (!pv->Within_sec_decl) && (!pv->Within_sos_decl) && (!pv->Within_free_decl)) {
    pv->Within_sec_decl = TRUE;
    pv->Within_sos_decl1 = FALSE;
  }
  return(SEC_SEC);
}

<INITIAL>{SOS} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  if(!pv->Within_sos_decl)
    pv->SOStype0 = (short)atoi(((char *)yytext) + 3);
  if((!pv->Within_int_decl) && (!pv->Within_sec_decl) && (!pv->Within_sos_decl) && (!pv->Within_free_decl))
    pv->Within_sos_decl = TRUE;
  return(SEC_SOS);
}

<INITIAL>{SOSD} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  FREE(pv->Last_var);
  pv->Last_var = strdup((char *)yytext);
  pv->Last_var[strlen(pv->Last_var) - 2] = 0;
  return(SOSDESCR);
}

<INITIAL>{FREE} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  if((!pv->Within_int_decl) && (!pv->Within_sec_decl) && (!pv->Within_sos_decl) && (!pv->Within_free_decl)) {
    pv->Within_free_decl = TRUE;
    pv->Within_sos_decl1 = FALSE;
  }
  return(SEC_FREE);
}

<INITIAL>{LBL} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;
  char *ptr;

  pp->lineno = yylineno;
  FREE(pv->Last_var);
  pv->Last_var = strdup((char *)yytext);
  ptr = pv->Last_var + strlen(pv->Last_var);
  ptr[-1] = ' ';
  while ((--ptr >= pv->Last_var) && (isspace(*ptr)))
    *ptr = 0;
  return(VARIABLECOLON);
}

<INITIAL>{VR} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  FREE(pv->Last_var);
  pv->Last_var = strdup((char *)yytext);
  return(VAR);
}

<INITIAL>":" {
  parse_parm *pp = PARM;

  pp->lineno = yylineno;
  return (COLON);
}

<INITIAL>"*" {
  parse_parm *pp = PARM;

  pp->lineno = yylineno;
  return(AR_M_OP);
}

<INITIAL>{OPEQ} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  pv->OP = *yytext;
  return(RE_OPEQ);
}

<INITIAL>{OPLE} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  pv->OP = *yytext;
  return(RE_OPLE);
}

<INITIAL>{OPGE} {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  pv->OP = *yytext;
  return(RE_OPGE);
}

<INITIAL>";" {
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pp->lineno = yylineno;
  pv->Within_int_decl = pv->Within_sec_decl = pv->Within_sos_decl = pv->Within_free_decl = FALSE;
  check_int_sec_sos_free_decl(pp, (int) pv->Within_int_decl, (int) pv->Within_sec_decl, (int) pv->Within_sos_decl, (int) pv->Within_free_decl);
  return(END_C);
}

<INITIAL>. {
  parse_parm *pp = PARM;

  pp->lineno = yylineno;
  report(NULL, CRITICAL, "LEX ERROR : %s lineno %d\n", yytext, yylineno);
  return(UNDEFINED);
}

%%
