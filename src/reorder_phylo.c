/* reorder_phylo.c       2008-03-17 */

/* Copyright 2008 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <R_ext/Applic.h>

void neworder_cladewise(int *n, int *edge1, int *edge2,
			int *N, int *neworder)
/* n: nb of tips, N: nb of edges */
{
    int i, j, k, node, *done, dn, *node_back, eb;
    done = &dn;
    node_back = &eb;

    /* done: indicates whether an edge has been collected
       node_back: the series of node from the root to `node'
       node: the current node */

    done = (int*)R_alloc(*N, sizeof(int));
    node_back = (int*)R_alloc(*N + 2 - *n, sizeof(int));
    memset(done, 0, *N * sizeof(int));

    j = k = 0;
    node = *n + 1;
    while (j < *N) {
        for (i = 0; i < *N; i++) {
  	    if (done[i] || edge1[i] != node) continue;
	    neworder[j] = i + 1;
	    j++;
	    done[i] = 1;
	    if (edge2[i] > *n) {
	        node_back[k] = node;
		k++;
		node = edge2[i];
		/* if found a new node, reset the loop */
		i = -1;
	    }
	}
	/* if arrived at the end of `edge', go down one node */
	k--;
	node = node_back[k];
    }
}

#define DO_NODE_PRUNING\
    /* go back down in `edge' to set `neworder' */\
    for (j = 0; j <= i; j++) {\
        /* if find the edge where `node' is */\
        /* the descendant, make as ready */\
        if (edge2[j] == node) ready[j] = 1;\
	if (edge1[j] != node) continue;\
	neworder[nextI] = j + 1;\
	ready[j] = 0; /* mark the edge as done */\
	nextI++;\
    }

void neworder_pruningwise(int *ntip, int *nnode, int *edge1,
			  int *edge2, int *nedge, int *neworder)
{
    int *Ndegr, degree, *ready, rdy, i, j, node, nextI, n;
    Ndegr = &degree;
    ready = &rdy;

    ready = (int*)R_alloc(*nedge, sizeof(int));

    /* use `nextI' temporarily because need an address for R_tabulate */
    nextI = *ntip +  *nnode;
    Ndegr = (int*)R_alloc(nextI, sizeof(int));
    memset(Ndegr, 0, nextI*sizeof(int));
    R_tabulate(edge1, nedge, &nextI, Ndegr);

    /* `ready' indicates whether an edge is ready to be */
    /* collected; only the terminal edges are initially ready */
    for (i = 0; i < *nedge; i++)
      if (edge2[i] <= *ntip) ready[i] = 1;
      else ready[i] = 0;

    /* `n' counts the number of times a node has been seen. */
    /* This algo will work if the tree is in cladewise order, */
    /* so that the nodes of "cherries" will be contiguous in `edge'. */
    n = 0;
    nextI = 0;
    while (nextI < *nedge - Ndegr[*ntip]) {
        for (i = 0; i < *nedge; i++) {
            if (!ready[i]) continue;
	    if (!n) {
	        /* if found an edge ready, initialize `node' and start counting */
	        node = edge1[i];
		n = 1;
	    } else { /* else counting has already started */
	        if (edge1[i] == node) n++;
		else {
		    /* if the node has changed we checked that all edges */
		    /* from `node' have been found */
		    if (n == Ndegr[node - 1]) {
		        DO_NODE_PRUNING
		    }
		    /* in all cases reset `n' and `node' and carry on */
		    node = edge1[i];
		    n = 1;
		}
	    } /* go to the next edge */
	    /* if at the end of `edge', check that we can't do a node */
	    if (n == Ndegr[node - 1]) {
	        DO_NODE_PRUNING
		n = 0;
	    }
        }
    }
    for (i = 0; i < *nedge; i++) {
        if (!ready[i]) continue;
	neworder[nextI] = i + 1;
	nextI++;
    }
}
