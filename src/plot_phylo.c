/* plot_phylo.c (2012-10-01) */

/* Copyright 2004-2012 Emmanuel Paradis

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>

void node_depth_edgelength(int *ntip, int *nnode, int *edge1, int *edge2,
			   int *nedge, double *edge_length, double *xx)
{
    int i;

    /* We do a preorder tree traversal starting from the bottom */
    /* of `edge'; we assume `xx' has 0 for the root and the tree */
    /* is in pruningwise order. */
    for (i = *nedge - 1; i >= 0; i--)
      xx[edge2[i] - 1] = xx[edge1[i] - 1] + edge_length[i];
}

void node_depth(int *ntip, int *nnode, int *edge1, int *edge2,
		int *nedge, double *xx)
{
    int i;

    /* First set the coordinates for all tips */
    for (i = 0; i < *ntip; i++) xx[i] = 1;

    /* Then compute recursively for the nodes; we assume `xx' has */
    /* been initialized with 0's which is true if it has been */
    /* created in R (the tree must be in pruningwise order) */
    for (i = 0; i < *nedge; i++)
      xx[edge1[i] - 1] = xx[edge1[i] - 1] + xx[edge2[i] - 1];
}

void node_height(int *ntip, int *nnode, int *edge1, int *edge2,
		 int *nedge, double *yy)
{
    int i, n;
    double S;

    /* The coordinates of the tips have been already computed */

    S = 0;
    n = 0;
    for (i = 0; i < *nedge; i++) {
	S += yy[edge2[i] - 1];
	n += 1;
        if (edge1[i + 1] != edge1[i]) {
	    yy[edge1[i] - 1] = S/n;
	    S = 0;
	    n = 0;
	}
    }
}

void node_height_clado(int *ntip, int *nnode, int *edge1, int *edge2,
		       int *nedge, double *xx, double *yy)
{
    int i, n;
    double S;

    node_depth(ntip, nnode, edge1, edge2, nedge, xx);

    /* The coordinates of the tips have been already computed */

    S = 0;
    n = 0;
    for (i = 0; i < *nedge; i++) {
	S += yy[edge2[i] - 1] * xx[edge2[i] - 1];
	n += xx[edge2[i] - 1];
        if (edge1[i + 1] != edge1[i]) {
	    yy[edge1[i] - 1] = S/n;
	    S = 0;
	    n = 0;
	}
    }
}
