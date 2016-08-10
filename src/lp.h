////////////////////////////////////////////////////////////////
//                                                            //
//  lp.h                                                      //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 21/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////


#ifndef __LP_h
#define __LP_h

#ifdef __cplusplus
extern "C" {
#endif

#include <gurobi_c.h>

#define CHECK_VAL_GRB(val,msg,env){                                                                 \
        if ((val)){                                                                                     \
            fprintf(stderr, "%s at %s, line %d: %s\n",(msg),__FILE__,__LINE__,GRBgeterrormsg (env) );   \
        }                                                                                               \
    }

#define CHECK_VAL_GRB2(val,msg,env){                                                                \
        if ((val)){                                                                                     \
            fprintf(stderr, "%s at %s, line %d: %s\n",(msg),__FILE__,__LINE__,GRBgeterrormsg (env) );   \
            goto CLEAN;                                                                                 \
        }                                                                                               \
        \
    }

typedef struct wctlp wctlp;

struct wctlp {
    GRBenv *env;
    GRBmodel *model;

    double dbl_cutoff;
};

typedef struct wctlp_warmstart {
    int rcount;
    int ccount;
    int *rstat;
    int *cstat;
    double *dnorm;
} wctlp_warmstart;



#define wctlp_CONT 0
#define wctlp_BIN 1
#define wctlp_INT 2

#define wctlp_EQUAL 'E'
#define wctlp_LESS_EQUAL 'L'
#define wctlp_GREATER_EQUAL 'G'

#define wctlp_LOWER 0
#define wctlp_BASIC 1
#define wctlp_UPPER 2
#define wctlp_FREE 3

#define wctlp_MIN 1
#define wctlp_MAX -1

int wctlp_init(wctlp **lp, const char *name);
void wctlp_free(wctlp **lp);

int wctlp_optimize(wctlp *lp, int *status);
int wctlp_objval(wctlp *, double *obj);
int wctlp_pi(wctlp *, double *pi);
int wctlp_x(wctlp *, double *x, int first);

int wctlp_basis_cols(wctlp *lp, int *cstat, int first);
int wctlp_change_obj(wctlp *lp, int start, int len, double *values);
int wctlp_addrow(wctlp *lp, int nzcount, int *cind, double *cval, char sense,
                 double rhs, char *name);
int wctlp_addcol(wctlp *lp, int nzcount, int *cind, double *cval, double obj,
                 double lb, double ub, char vartype, char *name);
int wctlp_deletecols(wctlp *lp, int first_cind, int last_cind);

int wctlp_set_coltypes(wctlp *lp, char sense);
int wctlp_setbound(wctlp *lp, int col, char lower_or_uper, double bound);
int wctlp_obj_sense(wctlp *lp, int sense);
int wctlp_setnodelimit(wctlp *lp, int mip_node_limit);
int wctlp_set_cutoff(wctlp *lp, double cutoff);

int wctlp_write(wctlp *lp, const char *fname);
void wctlp_warmstart_free(wctlp_warmstart **w);
void wctlp_printerrorcode(int c);

int wctlp_status(wctlp *lp, int *status);
int wctlp_chg_lb_var(wctlp *lp, int var, double lb);
int wctlp_pi_inf(wctlp *lp, double *pi);
int wctlp_get_nb_rows(wctlp *lp, int *nb_rows);

double lp_int_tolerance(void);

int wctlp_get_rhs(wctlp *lp, double *rhs);

#ifdef __cplusplus
}
#endif
#endif
