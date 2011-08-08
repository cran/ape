/* tree_build.c    2011-06-23 */

/* Copyright 2008-2011 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>

static int str2int(char *x, int n)
{
	int i, k = 1, ans = 0;

	for (i = n - 1; i >= 0; i--, k *= 10)
		ans += ((int)x[i] - 48) * k;

	return ans;
}

void extract_portion_Newick(const char *x, int a, int b, char *y)
{
	int i, j;

	for (i = a, j = 0; i <= b; i++, j++) y[j] = x[i];

	y[j] = '\0';
}

void decode_terminal_edge_token(const char *x, int a, int b, int *node, double *w)
{
	int co = a;
	char *endstr, str[100];

	while (x[co] != ':') co++;

	extract_portion_Newick(x, a, co - 1, str);
	*node = str2int(str, co - a);
	extract_portion_Newick(x, co + 1, b, str);
	*w = R_strtod(str, &endstr);
}

void decode_internal_edge(const char *x, int a, int b, char *lab, double *w)
{
	int co = a;
	char *endstr, str[100];

	while (x[co] != ':') co++;

	if (a == co) lab[0] = '\0'; /* if no node label */
	else extract_portion_Newick(x, a, co - 1, lab);

	extract_portion_Newick(x, co + 1, b, str);
	*w = R_strtod(str, &endstr);
}

#define ADD_INTERNAL_EDGE            \
    e[j] = curnode;                  \
    e[j + nedge] = curnode = ++node; \
    stack_internal[k++] = j;         \
    j++

#define ADD_TERMINAL_EDGE                                        \
    e[j] = curnode;                                              \
    decode_terminal_edge_token(x, pr + 1, ps - 1, &tmpi, &tmpd); \
    e[j + nedge] = tmpi;                                         \
    el[j] = tmpd;                                                \
    j++

#define GO_DOWN                                                  \
    decode_internal_edge(x, ps + 1, pt - 1, lab, &tmpd);         \
    SET_STRING_ELT(node_label, curnode - 1 - ntip, mkChar(lab)); \
    l = stack_internal[--k];					 \
    el[l] = tmpd;                                                \
    curnode = e[l]

SEXP treeBuildWithTokens(SEXP nwk)
{
	const char *x;
	int n, i, ntip = 1, nnode = 0, nedge, *e, curnode, node, j, *skeleton, nsk = 0, ps, pr, pt, tmpi, l, k, stack_internal[10000];
	double *el, tmpd;
	char lab[512];
	SEXP edge, edge_length, Nnode, node_label, phy;

	PROTECT(nwk = coerceVector(nwk, STRSXP));
	x = CHAR(STRING_ELT(nwk, 0));
	n = strlen(x);
	skeleton = (int *)R_alloc(n, sizeof(int *));

/* first pass on the Newick string to localize parentheses and commas */
	for (i = 0; i < n; i++) {
		if (x[i] == '(') {
			skeleton[nsk] = i;
			nsk++;
			continue;
		}
		if (x[i] == ',') {
			skeleton[nsk] = i;
			nsk++;
			ntip++;
			continue;
		}
		if (x[i] == ')') {
			skeleton[nsk] = i;
			nsk++;
			nnode++;
		}
	}

	nedge = ntip + nnode - 1;

	PROTECT(Nnode = allocVector(INTSXP, 1));
	PROTECT(edge = allocVector(INTSXP, nedge*2));
	PROTECT(edge_length = allocVector(REALSXP, nedge));
	PROTECT(node_label = allocVector(STRSXP, nnode));
	INTEGER(Nnode)[0] = nnode;

	e = INTEGER(edge);
	el = REAL(edge_length);

	curnode = node = ntip + 1;
	k = j = 0;
/* j: index of the current position in the edge matrix */
/* k: index of the current position in stack_internal */
/* stack_internal is a simple array storing the indices of the
   successive internal edges from the root; it's a stack so it is
   incremented every time an internal edge is added, and decremented
   every GO_DOWN step. This makes easy to the index of the
   subtending edge. */

/* second pass on the Newick string to build the "phylo" object elements */
	for (i = 1; i < nsk - 1; i++) {
		ps = skeleton[i];
//		Rprintf(""); /* <- again !!!! */
		if (x[ps] == '(') {
			ADD_INTERNAL_EDGE;
			continue;
		}
		pr = skeleton[i - 1];
		if (x[ps] == ',') {
			if (x[pr] != ')') {
				/* !!! accolades indispensables !!! */
				ADD_TERMINAL_EDGE;
			}
			continue;
		}
		if (x[ps] == ')') {
			pt = skeleton[i + 1]; // <- utile ???
			if (x[pr] == ',') {
				ADD_TERMINAL_EDGE;
 				GO_DOWN;
				continue;
			}
			if (x[pr] == ')') {
				GO_DOWN;
			}
		}
	}

	pr = skeleton[nsk - 2];
	ps = skeleton[nsk - 1];
/* is the last edge terminal? */
	if (x[pr] == ',' && x[ps] == ')') {
		ADD_TERMINAL_EDGE;
	}

/* is there a root edge and/or root label? */
	if (ps < n - 2) {
		i = ps + 1;
		while (i < n - 2 && x[i] != ':') i++;
		if (i < n - 2) {
			PROTECT(phy = allocVector(VECSXP, 5));
			SEXP root_edge;
			decode_internal_edge(x, ps + 1, n - 2, lab, &tmpd);
			PROTECT(root_edge = allocVector(REALSXP, 1));
			REAL(root_edge)[0] = tmpd;
			SET_VECTOR_ELT(phy, 4, root_edge);
			UNPROTECT(1);
			SET_STRING_ELT(node_label, 0, mkChar(lab));
		} else {
			extract_portion_Newick(x, ps + 1, n - 2, lab);
			SET_STRING_ELT(node_label, 0, mkChar(lab));
			PROTECT(phy = allocVector(VECSXP, 4));
		}
	} else PROTECT(phy = allocVector(VECSXP, 4));

	SET_VECTOR_ELT(phy, 0, edge);
	SET_VECTOR_ELT(phy, 1, edge_length);
	SET_VECTOR_ELT(phy, 2, Nnode);
	SET_VECTOR_ELT(phy, 3, node_label);

	UNPROTECT(6);
	return phy;
}

#undef ADD_INTERNAL_EDGE
#undef ADD_TERMINAL_EDGE
#undef GO_DOWN
