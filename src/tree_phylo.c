/* tree_phylo.c    2012-04-24 */

/* Copyright 2008-2012 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "me.h"

/* from newick.c */
int leaf(node *v);

static int curnod, curtip, iedge;

#define DO_EDGE\
	el[iedge] = EDGE->distance;\
	if (leaf(EDGE->head)) {\
		edge2[iedge] = curtip;\
		strncpy(tl[curtip - 1], EDGE->head->label, 6);\
		tl[curtip - 1][6] = '\0';\
		iedge++;\
		curtip++;\
	} else {\
		edge2[iedge] = curnod;\
		iedge++;\
		subtree2phylo(EDGE->head, edge1, edge2, el, tl);\
	}

void subtree2phylo(node *parent, int *edge1, int *edge2, double *el, char **tl)
{
	edge *EDGE; int localnode;

	EDGE = parent->leftEdge;
/* 'localnode' keeps a copy of the node ancestor # between
   the two (recursive) calls of subtree2phylo */
	localnode = edge1[iedge] = curnod;
	curnod++;
	DO_EDGE

	EDGE = parent->rightEdge;
	edge1[iedge] = localnode;
	DO_EDGE
}

/*
transforms a 'tree' struc of pointers into an object of class "phylo"
assumes the tree is unrooted and binary, so there are 2n - 3 edges
assumes labels are 6-char long
*/
void tree2phylo(tree *T, int *edge1, int *edge2, double *el, char **tl, int n)
{
	edge *EDGE;
	curnod = n + 1; /* the root for ape */

/* there's in fact only one edge from the "root" which is
   a tip in ape's terminology (i.e., a node of degree 1) */

	EDGE = T->root->leftEdge;
	edge1[0] = curnod;
	edge2[0] = 1; /* <- the 1st tip */
	strncpy(tl[0], T->root->label, 6);
	tl[0][6] = '\0';
	el[0] = EDGE->distance;
	/* now can initialize these two: */
	curtip = 2; /* <- the 2nd tip */
	iedge = 1; /* <- the 2nd edge */
	edge1[iedge] = curnod;

/* 'T->root->leftEdge->head' is the root for ape,
   so don't need to test if it's a leaf */

	subtree2phylo(EDGE->head, edge1, edge2, el, tl);
}
