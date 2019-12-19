/* ========================================================================= */
/* NAME  : lp_rlp.y                                                          */
/* ========================================================================= */

/*
   made reentrant with help of
   http://www.usualcoding.eu/post/2007/09/03/Building-a-reentrant-parser-in-C-with-Flex/Bison
*/

/*
   Note that a minimum version of bison is needed to be able to compile this.
   Older version don't know the reentrant code.
   Version 1.35 is not enough. v1.875 could be ok. Tested with v2.3
*/

%pure-parser
%parse-param {parse_parm *parm}
%parse-param {void *scanner}
%lex-param {yyscan_t *scanner}

%token VAR CONS INTCONS VARIABLECOLON INF SEC_INT SEC_BIN SEC_SEC SEC_SOS SOSDESCR SEC_FREE TOK_SIGN AR_M_OP RE_OPEQ RE_OPLE RE_OPGE END_C COMMA COLON MINIMISE MAXIMISE UNDEFINED


%{
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define scanner yyscanner
#define PARM yyget_extra(yyscanner)
#define YYSTYPE int
#define YY_EXTRA_TYPE parse_parm *
#define YY_FATAL_ERROR(msg) lex_fatal_error(PARM, yyscanner, msg)
#undef YY_INPUT
#define YY_INPUT(buf,result,max_size) result = lp_input((void *) PARM, buf, max_size);
#define yyerror read_error

#include "lpkit.h"
#include "yacc_read.h"

typedef struct parse_vars_s
{
  read_modeldata_func *lp_input;
  void *userhandle;
  char HadVar, HadVar0, HadVar1, HadVar2, HasAR_M_OP, HadConstraint, Had_lineair_sum, Had_lineair_sum0, do_add_row, HadSign, OP, Sign, isign, isign0, make_neg;
  char state, state0;
  char Within_int_decl;  /* TRUE when we are within an char declaration */
  char Within_bin_decl;  /* TRUE when we are within an bin declaration */
  char Within_sec_decl;  /* TRUE when we are within a sec declaration */
  char Within_sos_decl;  /* TRUE when we are within a sos declaration */
  char Within_sos_decl1;
  char Within_free_decl; /* TRUE when we are within a free declaration */
  short SOStype, SOStype0;        /* SOS type */
  int SOSNr;
  int SOSweight;         /* SOS weight */
  char *Last_var, *Last_var0;
  REAL f, f0, f1;
} parse_vars;

#ifdef FORTIFY
# include "lp_fortify.h"
#endif

/* let's please C++ users */
#ifdef __cplusplus
extern "C" {
#endif

#if defined MSDOS || defined __MSDOS__ || defined WINDOWS || defined _WINDOWS || defined WIN32 || defined _WIN32
#define YY_NO_UNISTD_H

static int isatty(int f)
{
  return(FALSE);
}

#if !defined _STDLIB_H
# define _STDLIB_H
#endif
#endif

static int __WINAPI lp_input_yyin(void *fpin, char *buf, int max_size)
{
  int result;

  result = fread( (char*)buf, sizeof(char), max_size, (FILE *)fpin);

  return(result);
}

static int __WINAPI lp_input(void *vpp, char *buf, int max_size)
{
  parse_parm *pp = (parse_parm *) vpp;
  parse_vars *pv = (parse_vars *) pp->parse_vars;
  int result;

  result = pv->lp_input(pv->userhandle, buf, max_size);
  if (result < 0)
    lex_fatal_error(pp, pp->scanner, "read() in flex scanner failed");
  return(result);
}

#ifdef __cplusplus
};
#endif

#include "lp_rlp.h"

#undef yylval

%}

%start inputfile
%%

EMPTY: /* EMPTY */
                ;

inputfile       :
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->isign = 0;
  pv->make_neg = 0;
  pv->Sign = 0;
  pv->HadConstraint = FALSE;
  pv->HadVar = pv->HadVar0 = FALSE;
}
                  objective_function
                  constraints
                  int_bin_sec_sos_free_declarations
                ;

/* start objective_function */

/*

 objective_function: MAXIMISE real_of | MINIMISE real_of | real_of;
 real_of:            lineair_sum END_C;
 lineair_sum:        EMPTY | x_lineair_sum;

*/

objective_function:   MAXIMISE real_of
{
  set_obj_dir(PARM, TRUE);
}
                    | MINIMISE real_of
{
  set_obj_dir(PARM, FALSE);
}
                    | real_of
                ;

real_of:            lineair_sum
                    END_C
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  add_row(pp);
  pv->HadConstraint = FALSE;
  pv->HadVar = pv->HadVar0 = FALSE;
  pv->isign = 0;
  pv->make_neg = 0;
}
                ;

lineair_sum:          EMPTY
                    | x_lineair_sum
                ;

/* end objective_function */



/* start constraints */

/*

 constraints:        EMPTY | x_constraints;
 x_constraints:      constraint | x_constraints constraint;
 constraint:         real_constraint | VARIABLECOLON real_constraint;
 real_constraint:    x_lineair_sum2 RE_OP x_lineair_sum3 optionalrange END_C;
 optionalrange:      EMPTY | RE_OP cons_term RHS_STORE;
 RE_OP:              RE_OPEQ | RE_OPLE | RE_OPGE;
 cons_term:          x_SIGN REALCONS | INF;
 x_lineair_sum2:     EMPTY | x_lineair_sum3;
 x_lineair_sum3:     x_lineair_sum | INF RHS_STORE;
 x_lineair_sum:      x_lineair_sum1;
 x_lineair_sum1:     x_lineair_term | x_lineair_sum1 x_lineair_term;
 x_lineair_term:     x_SIGN x_lineair_term1;
 x_lineair_term1:    REALCONS | optional_AR_M_OP VAR;
 x_SIGN:             EMPTY | TOK_SIGN;
 REALCONS:           INTCONS | CONS;
 optional_AR_M_OP:   EMPTY | AR_M_OP;

*/

constraints:      EMPTY
                | x_constraints
                ;

x_constraints   : constraint
                | x_constraints
                  constraint
                ;

constraint      : real_constraint
                | VARIABLECOLON
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if(!add_constraint_name(pp, pv->Last_var))
    YYABORT;
  pv->HadConstraint = TRUE;
}
                  real_constraint
                ;

real_constraint : x_lineair_sum2
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->HadVar1 = pv->HadVar0;
  pv->HadVar0 = FALSE;
}
                  RE_OP
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if(!store_re_op(pp, pv->OP, (int) pv->HadConstraint, (int) pv->HadVar, (int) pv->Had_lineair_sum))
    YYABORT;
  pv->make_neg = 1;
  pv->f1 = 0;
}
                  x_lineair_sum3
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->Had_lineair_sum0 = pv->Had_lineair_sum;
  pv->Had_lineair_sum = TRUE;
  pv->HadVar2 = pv->HadVar0;
  pv->HadVar0 = FALSE;
  pv->do_add_row = FALSE;
  if(pv->HadConstraint && !pv->HadVar ) {
    /* it is a range */
    /* already handled */
  }
  else if(!pv->HadConstraint && pv->HadVar) {
    /* it is a bound */

    if(!store_bounds(pp, TRUE))
      YYABORT;
  }
  else {
    /* it is a row restriction */
    if(pv->HadConstraint && pv->HadVar)
      store_re_op(pp, '\0', (int) pv->HadConstraint, (int) pv->HadVar, (int) pv->Had_lineair_sum); /* makes sure that data stored in temporary buffers is treated correctly */
    pv->do_add_row = TRUE;
  }
}
                  optionalrange
                  END_C
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if((!pv->HadVar) && (!pv->HadConstraint)) {
    yyerror(pp, pp->scanner, "parse error");
    YYABORT;
  }
  if(pv->do_add_row)
    add_row(pp);
  pv->HadConstraint = FALSE;
  pv->HadVar = pv->HadVar0 = FALSE;
  pv->isign = 0;
  pv->make_neg = 0;
  null_tmp_store(pp, TRUE);
}
                ;

optionalrange:    EMPTY
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if((!pv->HadVar1) && (pv->Had_lineair_sum0))
    if(!negate_constraint(pp))
      YYABORT;
}
                | RE_OP
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->make_neg = 0;
  pv->isign = 0;
  if(pv->HadConstraint)
    pv->HadVar = pv->Had_lineair_sum = FALSE;
  pv->HadVar0 = FALSE;
  if(!store_re_op(pp, (char) ((pv->OP == '<') ? '>' : (pv->OP == '>') ? '<' : pv->OP), (int) pv->HadConstraint, (int) pv->HadVar, (int) pv->Had_lineair_sum))
    YYABORT;
}
                  cons_term
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->f -= pv->f1;
}
                  RHS_STORE
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if((pv->HadVar1) || (!pv->HadVar2) || (pv->HadVar0)) {
    yyerror(pp, pp->scanner, "parse error");
    YYABORT;
  }

  if(pv->HadConstraint && !pv->HadVar ) {
    /* it is a range */
    /* already handled */
    if(!negate_constraint(pp))
      YYABORT;
  }
  else if(!pv->HadConstraint && pv->HadVar) {
    /* it is a bound */

    if(!store_bounds(pp, TRUE))
      YYABORT;
  }
}
                ;

x_lineair_sum2:   EMPTY
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  /* to allow a range */
  /* constraint: < max */
  if(!pv->HadConstraint) {
    yyerror(pp, pp->scanner, "parse error");
    YYABORT;
  }
  pv->Had_lineair_sum = FALSE;
}
                | x_lineair_sum3
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->Had_lineair_sum = TRUE;
}
                ;

x_lineair_sum3  :  x_lineair_sum
                | INF
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->isign = pv->Sign;
}
                  RHS_STORE
                ;

x_lineair_sum:
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->state = pv->state0 = 0;
}
                x_lineair_sum1
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if (pv->state == 1) {
    /* RHS_STORE */
    if (    (pv->isign0 || !pv->make_neg)
        && !(pv->isign0 && !pv->make_neg)) /* but not both! */
      pv->f0 = -pv->f0;
    if(pv->make_neg)
      pv->f1 += pv->f0;
    if(!rhs_store(pp, pv->f0, (int) pv->HadConstraint, (int) pv->HadVar, (int) pv->Had_lineair_sum))
      YYABORT;
  }
}
                ;

x_lineair_sum1  : x_lineair_term
                | x_lineair_sum1
                  x_lineair_term
                ;

x_lineair_term  : x_SIGN
                  x_lineair_term1
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if ((pv->HadSign || pv->state == 1) && (pv->state0 == 1)) {
    /* RHS_STORE */
    if (    (pv->isign0 || !pv->make_neg)
        && !(pv->isign0 && !pv->make_neg)) /* but not both! */
      pv->f0 = -pv->f0;
    if(pv->make_neg)
      pv->f1 += pv->f0;
    if(!rhs_store(pp, pv->f0, (int) pv->HadConstraint, (int) pv->HadVar, (int) pv->Had_lineair_sum))
      YYABORT;
  }
  if (pv->state == 1) {
    pv->f0 = pv->f;
    pv->isign0 = pv->isign;
  }
  if (pv->state == 2) {
    if((pv->HadSign) || (pv->state0 != 1)) {
     pv->isign0 = pv->isign;
     pv->f0 = 1.0;
    }
    if (    (pv->isign0 || pv->make_neg)
        && !(pv->isign0 && pv->make_neg)) /* but not both! */
      pv->f0 = -pv->f0;
    if(!var_store(pp, pv->Last_var, pv->f0, (int) pv->HadConstraint, (int) pv->HadVar, (int) pv->Had_lineair_sum)) {
      yyerror(pp, pp->scanner, "var_store failed");
      YYABORT;
    }
    pv->HadConstraint |= pv->HadVar;
    pv->HadVar = pv->HadVar0 = TRUE;
  }
  pv->state0 = pv->state;
}
                ;

x_lineair_term1 : REALCONS
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->state = 1;
}
                | optional_AR_M_OP
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if ((pv->HasAR_M_OP) && (pv->state != 1)) {
    yyerror(pp, pp->scanner, "parse error");
    YYABORT;
  }
}
                  VAR
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->state = 2;
}
                ;

RE_OP: RE_OPEQ | RE_OPLE | RE_OPGE
                ;

cons_term:        x_SIGN
                  REALCONS
                | INF
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->isign = pv->Sign;
}
                ;

/* end constraints */


/* start common for objective & constraints */

REALCONS: INTCONS | CONS
                ;

x_SIGN:           EMPTY
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->isign = 0;
  pv->HadSign = FALSE;
}
                | TOK_SIGN
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->isign = pv->Sign;
  pv->HadSign = TRUE;
}
                ;

optional_AR_M_OP: EMPTY
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->HasAR_M_OP = FALSE;
}
                | AR_M_OP
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->HasAR_M_OP = TRUE;
}
                ;

RHS_STORE:        EMPTY
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if (    (pv->isign || !pv->make_neg)
      && !(pv->isign && !pv->make_neg)) /* but not both! */
    pv->f = -pv->f;
  if(!rhs_store(pp, pv->f, (int) pv->HadConstraint, (int) pv->HadVar, (int) pv->Had_lineair_sum))
    YYABORT;
  pv->isign = 0;
}
                ;

/* end common for objective & constraints */



/* start int_bin_sec_sos_free_declarations */

int_bin_sec_sos_free_declarations:
                  EMPTY
                | real_int_bin_sec_sos_free_decls
                ;

real_int_bin_sec_sos_free_decls: int_bin_sec_sos_free_declaration
                | real_int_bin_sec_sos_free_decls int_bin_sec_sos_free_declaration
                ;

SEC_INT_BIN_SEC_SOS_FREE: SEC_INT | SEC_BIN | SEC_SEC | SEC_SOS | SEC_FREE
                ;

int_bin_sec_sos_free_declaration:
                  SEC_INT_BIN_SEC_SOS_FREE
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  pv->Within_sos_decl1 = pv->Within_sos_decl;
}
                  x_int_bin_sec_sos_free_declaration
                ;

xx_int_bin_sec_sos_free_declaration:
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if((!pv->Within_int_decl) && (!pv->Within_sec_decl) && (!pv->Within_sos_decl1) && (!pv->Within_free_decl)) {
    yyerror(pp, pp->scanner, "parse error");
    YYABORT;
  }
  pv->SOStype = pv->SOStype0;
  check_int_sec_sos_free_decl(pp, (int) pv->Within_int_decl, (int) pv->Within_sec_decl, (int) (pv->Within_sos_decl1 = (pv->Within_sos_decl1 ? 1 : 0)), (int) pv->Within_free_decl);
}
                  optionalsos
                  vars
                  optionalsostype
                  END_C
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if((pv->Within_sos_decl1) && (pv->SOStype == 0))
  {
    yyerror(pp, pp->scanner, "Unsupported SOS type (0)");
    YYABORT;
  }
}
                ;

x_int_bin_sec_sos_free_declaration:
                  xx_int_bin_sec_sos_free_declaration
                | x_int_bin_sec_sos_free_declaration xx_int_bin_sec_sos_free_declaration
                ;

optionalsos:      EMPTY
                | SOSDESCR
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  FREE(pv->Last_var0);
  pv->Last_var0 = strdup(pv->Last_var);
}
                  sosdescr
                ;

optionalsostype:  EMPTY
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if(pv->Within_sos_decl1) {
    set_sos_type(pp, pv->SOStype);
    set_sos_weight(pp, (double) pv->SOSweight, 1);
  }
}
                | RE_OPLE
                  INTCONS
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if((pv->Within_sos_decl1) && (!pv->SOStype))
  {
    set_sos_type(pp, pv->SOStype = (short) (pv->f + .1));
  }
  else
  {
    yyerror(pp, pp->scanner, "SOS type not expected");
    YYABORT;
  }
}
                optionalSOSweight
                ;

optionalSOSweight:EMPTY
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  set_sos_weight(pp, (double) pv->SOSweight, 1);
}
                | COLON
                  INTCONS
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  set_sos_weight(pp, pv->f, 1);
}
                ;

vars:             EMPTY
                | x_vars
                ;

x_vars          : onevarwithoptionalweight
                | x_vars
                  optionalcomma
                  onevarwithoptionalweight
                ;

optionalcomma:    EMPTY
                | COMMA
                ;

variable:         EMPTY
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if(pv->Within_sos_decl1 == 1)
  {
    char buf[16];

    pv->SOSweight++;
    sprintf(buf, "SOS%d", pv->SOSweight);
    storevarandweight(pp, buf);

    check_int_sec_sos_free_decl(pp, (int) pv->Within_int_decl, (int) pv->Within_sec_decl, 2, (int) pv->Within_free_decl);
    pv->Within_sos_decl1 = 2;
    pv->SOSNr = 0;
  }

  storevarandweight(pp, pv->Last_var);

  if(pv->Within_sos_decl1 == 2)
  {
    pv->SOSNr++;
    set_sos_weight(pp, (double) pv->SOSNr, 2);
  }
}
                ;

variablecolon:
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if(!pv->Within_sos_decl1) {
    yyerror(pp, pp->scanner, "parse error");
    YYABORT;
  }
  if(pv->Within_sos_decl1 == 1) {
    FREE(pv->Last_var0);
    pv->Last_var0 = strdup(pv->Last_var);
  }
  if(pv->Within_sos_decl1 == 2)
  {
    storevarandweight(pp, pv->Last_var);
    pv->SOSNr++;
    set_sos_weight(pp, (double) pv->SOSNr, 2);
  }
}
                ;

sosweight:        EMPTY
{
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if(pv->Within_sos_decl1 == 1)
  {
    char buf[16];

    pv->SOSweight++;
    sprintf(buf, "SOS%d", pv->SOSweight);
    storevarandweight(pp, buf);

    check_int_sec_sos_free_decl(pp, (int) pv->Within_int_decl, (int) pv->Within_sec_decl, 2, (int) pv->Within_free_decl);
    pv->Within_sos_decl1 = 2;
    pv->SOSNr = 0;

    storevarandweight(pp, pv->Last_var0);
    pv->SOSNr++;
  }

  set_sos_weight(pp, pv->f, 2);
}
                ;

sosdescr:         EMPTY
{ /* SOS name */
  parse_parm *pp = PARM;
  parse_vars *pv = (parse_vars *) pp->parse_vars;

  if(pv->Within_sos_decl1 == 1)
  {
    parse_parm *pp = PARM;
    parse_vars *pv = (parse_vars *) pp->parse_vars;

    storevarandweight(pp, pv->Last_var0);
    set_sos_type(pp, pv->SOStype);
    check_int_sec_sos_free_decl(pp, (int) pv->Within_int_decl, (int) pv->Within_sec_decl, 2, (int) pv->Within_free_decl);
    pv->Within_sos_decl1 = 2;
    pv->SOSNr = 0;
    pv->SOSweight++;
  }
}
                ;

onevarwithoptionalweight:
                  VAR
                  variable
                | VARIABLECOLON
                  variablecolon
                  INTCONSorVARIABLE
                ;

INTCONSorVARIABLE:REALCONS /* INTCONS */
                  sosweight
                | sosdescr
                  x_onevarwithoptionalweight
                ;

x_onevarwithoptionalweight:
                  VAR
                  variable
                | VARIABLECOLON
                  variablecolon
                  REALCONS /* INTCONS */
                  sosweight
                ;

/* end int_bin_sec_sos_free_declarations */

%%

static void yy_delete_allocated_memory(parse_parm *pp)
{
  parse_vars *pv = (parse_vars *) pp->parse_vars;
  /* free memory allocated by flex. Otherwise some memory is not freed.
     This is a bit tricky. There is not much documentation about this, but a lot of
     reports of memory that keeps allocated */

  /* If you get errors on this function call, just comment it. This will only result
     in some memory that is not being freed. */

# if defined YY_CURRENT_BUFFER
    /* flex defines the macro YY_CURRENT_BUFFER, so you should only get here if lp_rlp.h is
       generated by flex */
    /* lex doesn't define this macro and thus should not come here, but lex doesn't has
       this memory leak also ...*/

#  if 0
    /* older versions of flex */
    yy_delete_buffer(YY_CURRENT_BUFFER); /* comment this line if you have problems with it */
    yy_init = 1; /* make sure that the next time memory is allocated again */
    yy_start = 0;
#  else
    /* As of version 2.5.9 Flex  */
    yylex_destroy(pp->scanner); /* comment this line if you have problems with it */
#  endif
# endif

  FREE(pv->Last_var);
  FREE(pv->Last_var0);
}

static int parse(parse_parm *pp)
{
  return(yyparse(pp, pp->scanner));
}

lprec *read_lp1(lprec *lp, void *userhandle, read_modeldata_func read_modeldata, int verbose, char *lp_name)
{
  parse_vars *pv;
  lprec *lp1 = NULL;

  CALLOC(pv, 1, parse_vars);
  if (pv != NULL) {
    parse_parm pp;

    memset(&pp, 0, sizeof(pp));
    pp.parse_vars = (void *) pv;

    yylex_init(&pp.scanner);
    yyset_extra(&pp, pp.scanner);

    yyset_in((FILE *) userhandle, pp.scanner);
    yyset_out(NULL, pp.scanner);
    pv->lp_input = read_modeldata;
    pv->userhandle = userhandle;
    lp1 = yacc_read(lp, verbose, lp_name, parse, &pp, yy_delete_allocated_memory);
    FREE(pv);
  }
  return(lp1);
}

lprec * __WINAPI read_lp(FILE *filename, int verbose, char *lp_name)
{
  return(read_lp1(NULL, filename, lp_input_yyin, verbose, lp_name));
}

lprec * __WINAPI read_lpex(void *userhandle, read_modeldata_func read_modeldata, int verbose, char *lp_name)
{
  return(read_lp1(NULL, userhandle, read_modeldata, verbose, lp_name));
}

lprec *read_LP1(lprec *lp, char *filename, int verbose, char *lp_name)
{
  FILE *fpin;

  if((fpin = fopen(filename, "r")) != NULL) {
    lp = read_lp1(lp, fpin, lp_input_yyin, verbose, lp_name);
    fclose(fpin);
  }
  else
    lp = NULL;
  return(lp);
}

lprec * __WINAPI read_LP(char *filename, int verbose, char *lp_name)
{
  return(read_LP1(NULL, filename, verbose, lp_name));
}

MYBOOL __WINAPI LP_readhandle(lprec **lp, FILE *filename, int verbose, char *lp_name)
{
  if(lp != NULL)
    *lp = read_lp1(*lp, filename, lp_input_yyin, verbose, lp_name);

  return((lp != NULL) && (*lp != NULL));
}
