/* ape.c    2013-08-13 */

/* Copyright 2011-2013 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "ape.h"

int give_index(int i, int j, int n)
{
	if (i > j) return(DINDEX(j, i));
	else return(DINDEX(i, j));
}

/* declare functions here to register them below */

void C_additive(double *dd, int* np, int* mp, double *ret);
void BaseProportion(unsigned char *x, int *n, double *BF);
void C_bionj(double *X, int *N, int *edge1, int *edge2, double *el);
void C_bionjs(double *D, int *N, int *edge1, int *edge2, double *edge_length, int* fsS);
void delta_plot(double *D, int *size, int *nbins, int *counts, double *deltabar);
void dist_dna(unsigned char *x, int *n, int *s, int *model, double *d,
	      double *BF, int *pairdel, int *variance, double *var,
	      int *gamma, double *alpha);
void dist_nodes(int *n, int *m, int *e1, int *e2, double *el, int *N, double *D);
void C_ewLasso(double *D, int *N, int *e1, int *e2);
void GlobalDeletionDNA(unsigned char *x, int *n, int *s, int *keep);
void mat_expo(double *P, int *nr);
void me_b(double *X, int *N, int *labels,
	  int *nni, int *spr, int *tbr,
	  int *edge1, int *edge2, double *el);
void me_o(double *X, int *N, int *labels, int *nni,
	  int *edge1, int *edge2, double *el);
void C_mvr(double *D, double* v,int *N, int *edge1, int *edge2, double *edge_length);
void C_mvrs(double *D, double* v, int *N, int *edge1, int *edge2, double *edge_length, int* fsS);
void neworder_phylo(int *n, int *e1, int *e2, int *N, int *neworder, int *order);
void neworder_pruningwise(int *ntip, int *nnode, int *edge1,
			  int *edge2, int *nedge, int *neworder);
void C_nj(double *D, int *N, int *edge1, int *edge2, double *edge_length);
void C_njs(double *D, int *N, int *edge1, int *edge2, double *edge_length, int *fsS);
void node_depth(int *ntip, int *nnode, int *e1, int *e2,
		int *nedge, double *xx, int *method);
void node_depth_edgelength(int *ntip, int *nnode, int *edge1, int *edge2,
			   int *nedge, double *edge_length, double *xx);
void node_height(int *ntip, int *nnode, int *edge1, int *edge2,
		 int *nedge, double *yy);
void node_height_clado(int *ntip, int *nnode, int *edge1, int *edge2,
		       int *nedge, double *xx, double *yy);
void C_pic(int *ntip, int *nnode, int *edge1, int *edge2,
	   double *edge_len, double *phe, double *contr,
	   double *var_contr, int *var, int *scaled);
void C_rTraitCont(int *model, int *Nedge, int *edge1, int *edge2, double *el,
		  double *sigma, double *alpha, double *theta, double *x);
void SegSites(unsigned char *x, int *n, int *s, int *seg);
void C_treePop(int* splits, double* w,int* ncolp,int* np, int* ed1, int* ed2, double* edLen);
void C_triangMtd(double* d, int* np, int* ed1,int* ed2, double* edLen);
void C_triangMtds(double* d, int* np, int* ed1,int* ed2, double* edLen);
void C_ultrametric(double *dd, int* np, int* mp, double *ret);
void C_where(unsigned char *x, unsigned char *pat, int *s, int *p,
	   int *ans, int *n);

SEXP bipartition(SEXP edge, SEXP nbtip, SEXP nbnode);
SEXP prop_part(SEXP TREES, SEXP nbtree, SEXP keep_partitions);
SEXP rawStreamToDNAbin(SEXP x);
SEXP seq_root2tip(SEXP edge, SEXP nbtip, SEXP nbnode);
SEXP treeBuildWithTokens(SEXP nwk);

static R_CMethodDef C_entries[] = {
    {"C_additive", (DL_FUNC) &C_additive, 4},
    {"BaseProportion", (DL_FUNC) &BaseProportion, 3},
    {"C_bionj", (DL_FUNC) &C_bionj, 5},
    {"C_bionjs", (DL_FUNC) &C_bionjs, 6},
    {"delta_plot", (DL_FUNC) &delta_plot, 5},
    {"dist_dna", (DL_FUNC) &dist_dna, 11},
    {"dist_nodes", (DL_FUNC) &dist_nodes, 7},
    {"C_ewLasso", (DL_FUNC) &C_ewLasso, 4},
    {"GlobalDeletionDNA", (DL_FUNC) &GlobalDeletionDNA, 4},
    {"mat_expo", (DL_FUNC) &mat_expo, 2},
    {"me_b", (DL_FUNC) &me_b, 9},
    {"me_o", (DL_FUNC) &me_o, 7},
    {"C_mvr", (DL_FUNC) &C_mvr, 6},
    {"C_mvrs", (DL_FUNC) &C_mvrs, 7},
    {"neworder_phylo", (DL_FUNC) &neworder_phylo, 6},
    {"neworder_pruningwise", (DL_FUNC) &neworder_pruningwise, 6},
    {"C_nj", (DL_FUNC) &C_nj, 5},
    {"C_njs", (DL_FUNC) &C_njs, 6},
    {"node_depth", (DL_FUNC) &node_depth, 7},
    {"node_depth_edgelength", (DL_FUNC) &node_depth_edgelength, 7},
    {"node_height", (DL_FUNC) &node_height, 6},
    {"node_height_clado", (DL_FUNC) &node_height_clado, 7},
    {"C_pic", (DL_FUNC) &C_pic, 10},
    {"C_rTraitCont", (DL_FUNC) &C_rTraitCont, 9},
    {"SegSites", (DL_FUNC) &SegSites, 4},
    {"C_treePop", (DL_FUNC) &C_treePop, 7},
    {"C_triangMtd", (DL_FUNC) &C_triangMtd, 5},
    {"C_triangMtds", (DL_FUNC) &C_triangMtds, 5},
    {"C_ultrametric", (DL_FUNC) &C_ultrametric, 4},
    {"C_where", (DL_FUNC) &C_where, 6},
    {NULL, NULL, 0}
};

/*
static R_CMethodDef C_entries[] = {
    {"BaseProportion", (DL_FUNC) &BaseProportion, 3, {RAWSXP, INTSXP, REALSXP}},
    {"C_bionj", (DL_FUNC) &C_bionj, 5, {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP}},
    {"delta_plot", (DL_FUNC) &delta_plot, 5, {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP}},
    {"dist_dna", (DL_FUNC) &dist_dna, 11, {RAWSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP}},
    {"dist_nodes", (DL_FUNC) &dist_nodes, 7, {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP}},
    {"C_ewLasso", (DL_FUNC) &C_ewLasso, 4, {REALSXP, INTSXP, INTSXP, INTSXP}},
    {"GlobalDeletionDNA", (DL_FUNC) &GlobalDeletionDNA, 4, {RAWSXP, INTSXP, INTSXP, INTSXP}},
    {"mat_expo", (DL_FUNC) &mat_expo, 2, {REALSXP, INTSXP}},
    {"me_b", (DL_FUNC) &me_b, 9, {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP}},
    {"me_o", (DL_FUNC) &me_o, 7, {REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP}},
    {"C_mvr", (DL_FUNC) &C_mvr, 6, {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP}},
    {"C_mvrs", (DL_FUNC) &C_mvrs, 7, {REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP}},
    {"neworder_phylo", (DL_FUNC) &neworder_phylo, 6, {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP}},
    {"neworder_pruningwise", (DL_FUNC) &neworder_pruningwise, 6, {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP}},
    {"C_nj", (DL_FUNC) &C_nj, 5, {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP}},
    {"node_depth", (DL_FUNC) &node_depth, 7, {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP}},
    {"node_depth_edgelength", (DL_FUNC) &node_depth_edgelength, 7, {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP}},
    {"node_height", (DL_FUNC) &node_height, 6, {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP}},
    {"node_height_clado", (DL_FUNC) &node_height_clado, 7, {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP}},
    {"C_pic", (DL_FUNC) &C_pic, 9, {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP}},
    {"C_rTraitCont", (DL_FUNC) &C_rTraitCont, 9, {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP}},
    {"SegSites", (DL_FUNC) &SegSites, 4, {RAWSXP, INTSXP, INTSXP, INTSXP}},
    {"C_treePop", (DL_FUNC) &C_treePop, 7, {INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP}},
    {"C_triangMtd", (DL_FUNC) &C_triangMtd, 5, {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP}},
    {"C_triangMtds", (DL_FUNC) &C_triangMtds, 5, {REALSXP, INTSXP, INTSXP, INTSXP, REALSXP}},
    {"C_where", (DL_FUNC) &C_where, 6, {RAWSXP, RAWSXP, INTSXP, INTSXP, INTSXP, INTSXP}},
    {NULL, NULL, 0}
};
*/

static R_CallMethodDef Call_entries[] = {
    {"bipartition", (DL_FUNC) &bipartition, 3},
    {"prop_part", (DL_FUNC) &prop_part, 3},
    {"rawStreamToDNAbin", (DL_FUNC) &rawStreamToDNAbin, 1},
    {"seq_root2tip", (DL_FUNC) &seq_root2tip, 3},
    {"treeBuildWithTokens", (DL_FUNC) &treeBuildWithTokens, 1},
    {NULL, NULL, 0}
};

void R_init_ape(DllInfo *info)
{
    R_registerRoutines(info, C_entries, Call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
