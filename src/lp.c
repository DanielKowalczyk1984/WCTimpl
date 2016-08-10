////////////////////////////////////////////////////////////////
//                                                            //
//  lp.c                                                      //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 21/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////


#include "lp.h"
#include "util.h"





const double int_tolerance = 0.00001;


int wctlp_init(wctlp **lp, const char *name) {
    int val = 0;
    (*lp) = (wctlp *) CC_SAFE_MALLOC(1, wctlp);
    CCcheck_NULL(lp, "Out of memory for lp");
    (*lp)->env = (GRBenv *) NULL;
    (*lp)->model = (GRBmodel *) NULL;
    val = GRBloadenv(&((*lp)->env), NULL);

    if (val) {
        printf("val = %d\n", val);
    }

    CCcheck_val(val, "GRBloadenv failed");
    val = GRBsetintparam((*lp)->env, GRB_INT_PAR_OUTPUTFLAG,
                         (dbg_lvl() > 1) ? 1 : 0);
    CHECK_VAL_GRB(val, "GRBsetintparam OUTPUTFLAG failed", (*lp)->env);
    val = GRBsetintparam((*lp)->env, GRB_INT_PAR_THREADS, 1);
    CHECK_VAL_GRB(val, "GRBsetintparam TREADS failed", (*lp)->env);
    val = GRBsetintparam((*lp)->env, GRB_INT_PAR_METHOD, GRB_METHOD_PRIMAL);
    CHECK_VAL_GRB(val, "GRBsetintparam LPMETHOD failed", (*lp)->env);
    val = GRBsetintparam((*lp)->env, GRB_INT_PAR_INFUNBDINFO, 1);
    CHECK_VAL_GRB(val, "GRBsetintparam INFUNBDINFO", (*lp)->env);
    val =  GRBsetintparam((*lp)->env, GRB_INT_PAR_PRESOLVE, 0);
    CHECK_VAL_GRB(val, "Failed in setting parameter", (*lp)->env);
    val = GRBnewmodel((*lp)->env, &((*lp)->model), name, 0, (double *)NULL,
                      (double *)NULL, (double *)NULL, (char *)NULL, NULL);
    CHECK_VAL_GRB(val, "GRBnewmodel failed", (*lp)->env);
    return val;
}

void wctlp_free(wctlp **lp) {
    if (*lp) {
        if ((*lp)->model) {
            GRBfreemodel((*lp)->model);
        }

        if ((*lp)->env) {
            GRBfreeenv((*lp)->env);
        }

        free(*lp);
        *lp = (wctlp *)NULL;
    }
}

int wctlp_optimize(wctlp *lp, int *status) {
    int val;
    val = GRBoptimize((lp)->model);
    CHECK_VAL_GRB(val, "GRBoptimize failed", lp->env);
    val = GRBgetintattr(lp->model, GRB_INT_ATTR_STATUS, status);
    CHECK_VAL_GRB(val, "GRBgetinattr failed", lp->env);

    if (*status != GRB_OPTIMAL && *status != GRB_INFEASIBLE) {
        printf("Failed to solve the model to optimality. status= ");

        switch (*status) {
        case GRB_LOADED:
            printf("GRB_LOADED");
            val = 1;
            break;

        case GRB_UNBOUNDED:
            printf("GRB_UNBOUNDED");
            val = 1;
            break;

        case GRB_INF_OR_UNBD:
            printf("GRB_INF_OR_UNBD");
            val = 1;
            break;

        default:
            printf("%d", *status);
        }

        printf("\n");

        if (val) {
            goto CLEAN;
        }
    }

CLEAN:
    return val;
}

int wctlp_status(wctlp *lp, int *status) {
    int val = 0;
    val = GRBgetintattr(lp->model, GRB_INT_ATTR_STATUS, status);
    CHECK_VAL_GRB(val, "GRBgetintattr failed", lp->env);
    return val;
}

int wctlp_objval(wctlp *lp, double *obj) {
    int val = 0;
    val = GRBgetdblattr(lp->model, GRB_DBL_ATTR_OBJVAL, obj);
    CHECK_VAL_GRB(val, "GRBgetdblattr OBJVAL failed", lp->env);
    return val;
}

int wctlp_change_obj(wctlp *lp, int start, int len, double *values) {
    int val = 0;
    val = GRBsetdblattrarray(lp->model, GRB_DBL_ATTR_OBJ, start, len, values);
    CHECK_VAL_GRB(val, "Failed in GRBsetdblattrarray", lp->env);
    val = GRBupdatemodel(lp->model);
    CHECK_VAL_GRB(val, "Failed to update model", lp->env);
    return val;
}


int wctlp_addrow(wctlp *lp, int nzcount, int *cind , double *cval,
                 char sense, double rhs, char *name) {
    int val = 0;
    char isense;

    switch (sense) {
    case wctlp_EQUAL:
        isense = GRB_EQUAL;
        break;

    case wctlp_LESS_EQUAL:
        isense = GRB_LESS_EQUAL;
        break;

    case wctlp_GREATER_EQUAL:
        isense = GRB_GREATER_EQUAL;
        break;

    default:
        fprintf(stderr, "Unknown variable sense: %c\n", sense);
        val = 1;
        return val;
    }

    val = GRBaddconstr(lp->model, nzcount, cind, cval, isense, rhs, name);
    CHECK_VAL_GRB(val, "Failed GRBadd", lp->env);
    val = GRBupdatemodel(lp->model);
    CHECK_VAL_GRB(val, "Failed updating the model", lp->env);
    return val;
}
int wctlp_addcol(wctlp *lp, int nzcount, int *cind , double *cval,
                 double obj, double lb, double ub, char sense, char *name) {
    int val = 0;
    char isense;

    switch (sense) {
    case wctlp_CONT:
        isense = GRB_CONTINUOUS;
        break;

    case wctlp_BIN:
        isense = GRB_BINARY;
        break;

    case wctlp_INT:
        isense = GRB_INTEGER;
        break;

    default:
        fprintf(stderr, "Unknown variable sense: %c\n", sense);
        val = 1;
        return val;
    }

    val = GRBaddvar(lp->model, nzcount, cind, cval, obj, lb, ub, isense, name);
    CHECK_VAL_GRB(val, "Failed adding GRBaddvar", lp->env);
    val = GRBupdatemodel(lp->model);
    CHECK_VAL_GRB(val, "Failed updating the model", lp->env);
    return val;
}

int wctlp_deletecols(wctlp *lp, int first, int last) {
    int val = 0;
    int *dellist = (int *) NULL;
    int ndel = last - first + 1 ;
    int i;
    dellist = CC_SAFE_MALLOC(ndel, int);
    CCcheck_NULL_2(dellist, "Failed to allocated memory to dellist");

    for (i = 0; i < ndel; ++i) {
        dellist[i] = first + i;
    }

    val = GRBdelvars(lp->model, ndel, dellist);
    CHECK_VAL_GRB2(val, "Failed to delete cols", lp->env);
    GRBupdatemodel(lp->model);
    CHECK_VAL_GRB2(val, "Failed to update the model", lp->env);
CLEAN:
    CC_IFFREE(dellist, int);
    return val;
}

int wctlp_pi(wctlp *lp, double *pi) {
    int val = 0;
    int nrows;
    int solstat;
    val = GRBgetintattr(lp->model, GRB_INT_ATTR_STATUS, &solstat);
    CHECK_VAL_GRB(val, "failed to get attribute status gurobi", lp->env);
    val = GRBgetintattr(lp->model, GRB_INT_ATTR_NUMCONSTRS, &nrows);
    CHECK_VAL_GRB(val, "Failed to get nrows", lp->env);

    if (nrows == 0) {
        fprintf(stderr, "No rows in the LP\n");
        val = 1;
        return val;
    }

    val = GRBgetdblattrarray(lp->model, GRB_DBL_ATTR_PI, 0, nrows, pi);
    CHECK_VAL_GRB(val, "Failed to get the dual prices", lp->env);
    return val;
}

int wctlp_pi_inf(wctlp *lp, double *pi) {
    int val = 0;
    int nrows;
    int solstat;
    val = GRBgetintattr(lp->model, GRB_INT_ATTR_STATUS, &solstat);
    CHECK_VAL_GRB(val, "failed to get attribute status gurobi", lp->env);
    val = GRBgetintattr(lp->model, GRB_INT_ATTR_NUMCONSTRS, &nrows);
    CHECK_VAL_GRB(val, "Failed to get nrows", lp->env);

    if (nrows == 0) {
        fprintf(stderr, "No rows in the LP\n");
        val = 1;
        return val;
    }

    val = GRBgetdblattrarray(lp->model, GRB_DBL_ATTR_FARKASDUAL, 0, nrows, pi);
    CHECK_VAL_GRB(val, "Failed to get the dual prices", lp->env);
    return val;
}

int wctlp_x(wctlp *lp, double *x, int first) {
    int val = 0;
    int ncols;
    int solstat;
    val = GRBgetintattr(lp->model, GRB_INT_ATTR_STATUS, &solstat);
    CHECK_VAL_GRB(val, "Failed to the status of model", lp->env);

    if (solstat == GRB_INFEASIBLE) {
        fprintf(stderr, "Problem is infeasible\n");
        val = 1;
        return val;
    }

    val = GRBgetintattr(lp->model, GRB_INT_ATTR_NUMVARS , &ncols);
    CHECK_VAL_GRB(val, "", lp->env);

    if (ncols == 0) {
        fprintf(stderr, "Lp has no variables\n");
        val = 1;
        return val;
    }

    val = GRBgetdblattrarray(lp->model, GRB_DBL_ATTR_X, first, ncols - first, x);
    CHECK_VAL_GRB(val, "Failed in GRB_DBL_ATTR_X", lp->env);
    return val;
}

int wctlp_basis_cols(wctlp *lp, int *cstat, int first) {
    int val = 0;
    int ncols, i;
    val = GRBgetintattr(lp->model, GRB_INT_ATTR_NUMVARS, &ncols);
    CHECK_VAL_GRB(val, "Failed to get", lp->env);
    val = GRBgetintattrarray(lp->model, GRB_INT_ATTR_VBASIS, first, ncols - first,
                             cstat);
    CHECK_VAL_GRB(val, "Failed to get", lp->env);

    for (i = 0; i < ncols - first; ++i) {
        switch (cstat[i]) {
        case GRB_BASIC:
            cstat[i] = wctlp_BASIC;
            break;

        case GRB_NONBASIC_LOWER:
            cstat[i] = wctlp_LOWER;
            break;

        case GRB_NONBASIC_UPPER:
            cstat[i] = wctlp_UPPER;
            break;

        case GRB_SUPERBASIC:
            cstat[i] = wctlp_FREE;
            break;

        default:
            val = 1;
            CHECK_VAL_GRB(val, "ERROR: Received unknwn cstat", lp->env);
        }
    }

    return val;
}

int wctlp_set_coltypes(wctlp *lp, char sense) {
    int nvars, i, val = 0;
    char isense;

    switch (sense) {
    case wctlp_CONT:
        isense = GRB_CONTINUOUS;
        break;

    case wctlp_BIN:
        isense = GRB_BINARY;
        break;

    case wctlp_INT:
        isense = GRB_INTEGER;
        break;

    default:
        fprintf(stderr, "Unknown variable sense: %c\n", sense);
        val = 1;
        return val;
    }

    val = GRBgetintattr(lp->model, GRB_INT_ATTR_NUMVARS, &nvars);
    CHECK_VAL_GRB(val, "Failed to get number of variables", lp->env);

    for (i = 0; i < nvars; ++i) {
        val = GRBsetcharattrelement(lp->model, GRB_CHAR_ATTR_VTYPE, i, isense);
        CHECK_VAL_GRB(val , "Failed to set variable types", lp->env);
    }

    val = GRBupdatemodel(lp->model);
    CHECK_VAL_GRB(val, "Failed to update model", lp->env);
    return val;
}

int wctlp_get_rhs(wctlp *lp, double *rhs) {
    int val = 0;
    int nconstr;
    val = GRBgetintattr(lp->model, GRB_INT_ATTR_NUMCONSTRS, &nconstr);
    CHECK_VAL_GRB(val, "Failed in getting the number of variables", lp->env);
    val = GRBgetdblattrarray(lp->model, GRB_DBL_ATTR_RHS, 0, nconstr, rhs);
    CHECK_VAL_GRB(val, "Failed in getting the RHS", lp->env);
    return  val;
}


int wctlp_setbound(wctlp *lp, int col, char lb_or_ub, double bound) {
    int val = 0;

    if (lb_or_ub == 'L') {
        val = GRBsetdblattrelement(lp->model, GRB_DBL_ATTR_LB, col, bound);
    } else {
        val = GRBsetdblattrelement(lp->model, GRB_DBL_ATTR_UB, col, bound);
    }

    CHECK_VAL_GRB(val, "Failed to set bound", lp->env);
    val = GRBupdatemodel(lp->model);
    CHECK_VAL_GRB(val, "Failed to update model", lp->env);
    return val;
}

int wctlp_obj_sense(wctlp *lp, int sense) {
    int val = 0;
    val = GRBsetintattr(lp->model, GRB_INT_ATTR_MODELSENSE, sense);
    CHECK_VAL_GRB(val, "Failed to set obj sense", lp->env);
    return val;
}

int wctlp_write(wctlp *lp, const char *fname) {
    int val  = 0;
    val = GRBwrite(lp->model, fname);
    CHECK_VAL_GRB(val, "failed GRBwrite", lp->env);
    return val;
}

int wctlp_setnodelimit(wctlp *lp, int mip_node_limit) {
    int rval = GRBsetdblparam(GRBgetenv(lp->model), GRB_DBL_PAR_NODELIMIT,
                              mip_node_limit);
    CHECK_VAL_GRB2(rval, "GRBsetdblparam NODELIMIT failed", lp->env);
CLEAN:
    return rval;
}

void wctlp_warmstart_free(wctlp_warmstart **w) {
    if (*w != (wctlp_warmstart *) NULL) {
        CC_IFFREE((*w)-> cstat, int);
        CC_IFFREE((*w)-> rstat, int);
        CC_IFFREE((*w)->dnorm, double);
        CC_FREE(*w, wctlp_warmstart);
    }
}

int wctlp_get_nb_rows(wctlp *lp, int *nb_rows) {
    int val = 0;
    val = GRBgetintattr(lp->model, GRB_INT_ATTR_NUMCONSTRS, nb_rows);
    CHECK_VAL_GRB(val, "Failed to get the number of variables", lp->env);
    return val;
}



void wctlp_printerrorcode(int c) {
    switch (c) {
    case GRB_ERROR_OUT_OF_MEMORY:
        printf("Available memory was exhausted\n");
        break;

    case GRB_ERROR_NULL_ARGUMENT:
        printf("NULL input valua provided for a required argument\n");
        break;

    case GRB_ERROR_INVALID_ARGUMENT:
        printf("An invalid value was provided for a routine argument\n");
        break;

    case GRB_ERROR_UNKNOWN_ATTRIBUTE:
        printf("Tried to query or set an unknown attribute\n");
        break;

    case GRB_ERROR_DATA_NOT_AVAILABLE:
        printf("Attempted to query or set an attribute that could\n");
        printf("not be acessed at that time\n");
        break;

    case GRB_ERROR_UNKNOWN_PARAMETER:
        printf("Tried to query or set an unknown parameter\n");
        break;

    case GRB_ERROR_VALUE_OUT_OF_RANGE:
        printf("Tried to set a parameter to value that is outside\n");
        printf("the parameter's range\n");
        break;

    case GRB_ERROR_NO_LICENSE:
        printf("Failed to obtain a valid license\n");
        break;

    case GRB_ERROR_SIZE_LIMIT_EXCEEDED:
        printf("Attempted to solve a model that is larger than the\n");
        printf("limit for a demo version\n");
        break;

    case GRB_ERROR_CALLBACK:
        printf("Problem with callback\n");
        break;

    case GRB_ERROR_FILE_READ:
        printf("Failed to read the requested file\n");
        break;

    case GRB_ERROR_FILE_WRITE:
        printf("Failed to write the requested file\n");
        break;

    case GRB_ERROR_NUMERIC:
        printf("Numerical error during the requested operation\n");
        break;

    case GRB_ERROR_IIS_NOT_INFEASIBLE:
        printf("Attempted to perform infeasibility analysis on a feasible model\n");
        break;

    case GRB_ERROR_NOT_FOR_MIP:
        printf("Requested model not valid for a mip model\n");
        break;

    case GRB_ERROR_OPTIMIZATION_IN_PROGRESS:
        printf("Tried to query or modify a model while optimization was in progress\n");
        break;

    default:
        printf("Unknown error code: %d", c);
    }
}

double lp_int_tolerance() {
    return int_tolerance;
}
