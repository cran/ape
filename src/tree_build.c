/* tree_build.c    2009-11-21 */

/* Copyright 2008-2009 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>

static int str2int(char *x, int n)
{
	int i, k = 1, ans = 0;

	for (i = n - 2; i >= 0; i--, k *= 10)
		ans += ((int)x[i] - 48) * k;

	return ans;
}

void decode_edge(const char *x, int a, int b, int *node, double *w)
{
	int i, k, co = a;
	char *endstr, str[100];

	while (x[co] != ':') co++;
	if (a == co) *node = 0;
	else {
		for (i = a, k = 0; i < co; i++, k++) str[k] = x[i];
		str[k] = '\0';
		*node = str2int(str, k + 1);
	}
	for (i = co + 1, k = 0; i <= b; i++, k++) str[k] = x[i];
	str[k] = '\0';
	*w = R_strtod(str, &endstr);
}

#define ADD_INTERNAL_EDGE\
    e[j] = curnode;\
    e[j + nedge] = curnode = ++node;\
    j++

#define ADD_TERMINAL_EDGE\
    e[j] = curnode;\
    decode_edge(x, pr + 1, ps - 1, &tmpi, &tmpd);\
    e[j + nedge] = tmpi;\
    el[j] = tmpd;\
    j++

#define GO_DOWN\
    l = j - 1;\
    while (e[l + nedge] != curnode) l--;\
    decode_edge(x, ps + 1, pt - 1, &tmpi, &tmpd);\
    el[l] = tmpd;\
    curnode = e[l]

SEXP treeBuildWithTokens(SEXP nwk)
{
	const char *x;
	int n, i, ntip = 1, nnode = 0, nedge, *e, curnode, node, j, *skeleton, nsk = 0, ps, pr, pt, tmpi, l;
	double *el, tmpd;
	SEXP edge, edge_length, Nnode, phy;

	PROTECT(nwk = coerceVector(nwk, STRSXP));
	x = CHAR(STRING_ELT(nwk, 0));
	n = strlen(x);
	skeleton = (int *)R_alloc(n, sizeof(int *));
	for (i = 0; i < n; i++) {
		if (x[i] == '(' || x[i] == ',' || x[i] == ')') {
			skeleton[nsk] = i;
			nsk++;
			switch(x[i]) {
			case '(': break;
			case ',': ntip++; break;
			case ')': nnode++; break;
			}
		}
	}
	nedge = ntip + nnode - 1;

	PROTECT(Nnode = allocVector(INTSXP, 1));
	PROTECT(edge = allocVector(INTSXP, nedge*2));
	PROTECT(edge_length = allocVector(REALSXP, nedge));
	INTEGER(Nnode)[0] = nnode;

	e = INTEGER(edge);
	el = REAL(edge_length);

	curnode = node = ntip + 1;
	j = 0;

	for (i = 1; i < nsk - 1; i++) {
		ps = skeleton[i];
		Rprintf(""); /* <- again !!!! */
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
			pt = skeleton[i + 1];
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

/* is there a root edge? */
	if (ps < n - 2) {
		PROTECT(phy = allocVector(VECSXP, 4));
		SEXP root_edge;
		decode_edge(x, ps + 1, n - 2, &tmpi, &tmpd);
		PROTECT(root_edge = allocVector(REALSXP, 1));
		REAL(root_edge)[0] = tmpd;
		SET_VECTOR_ELT(phy, 3, root_edge);
		UNPROTECT(1);
	} else PROTECT(phy = allocVector(VECSXP, 3));

	SET_VECTOR_ELT(phy, 0, edge);
	SET_VECTOR_ELT(phy, 1, edge_length);
	SET_VECTOR_ELT(phy, 2, Nnode);

	UNPROTECT(5);
	return phy;
}

#undef ADD_INTERNAL_EDGE
#undef ADD_TERMINAL_EDGE
#undef GO_DOWN
