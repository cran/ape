/* bipartition.c    2006-10-06 */

/* Copyright 2005-2006 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>

SEXP seq_root2tip(SEXP edge1, SEXP edge2, SEXP nbtip, SEXP nbnode)
{
    int i, j, k, nbedge, *x1, *x2, *done, dn, sumdone, lt, ROOT, Ntip, Nnode;
    SEXP ans, seqnod, tmp_vec;

    PROTECT(edge1 = coerceVector(edge1, INTSXP));
    PROTECT(edge2 = coerceVector(edge2, INTSXP));
    PROTECT(nbtip = coerceVector(nbtip, INTSXP));
    PROTECT(nbnode = coerceVector(nbnode, INTSXP));
    x1 = INTEGER(edge1);
    x2 = INTEGER(edge2);
    Ntip = *INTEGER(nbtip);
    Nnode = *INTEGER(nbnode);
    nbedge = LENGTH(edge1);
    ROOT = Ntip + 1;

    PROTECT(ans = allocVector(VECSXP, Ntip));
    PROTECT(seqnod = allocVector(VECSXP, Nnode));

    done = &dn;
    done = (int*)R_alloc(Nnode, sizeof(int));
    for (i = 0; i < Nnode; i++) done[i] = 0;

    tmp_vec = allocVector(INTSXP, 1);
    INTEGER(tmp_vec)[0] = ROOT; /* sure ? */
    SET_VECTOR_ELT(seqnod, 0, tmp_vec);
    sumdone = 0;

    while (sumdone < Nnode) {
        for (i = 0; i < Nnode; i++) { /* loop through all nodes */
	    /* if the vector is not empty and its */
	    /* descendants are not yet found */
	    if (VECTOR_ELT(seqnod, i) == R_NilValue || done[i]) continue;
	    /* look for the descendants in 'edge': */
	    for (j = 0; j < nbedge; j++) {
	        /* skip the terminal edges, we look only for nodes now */
	        if (x1[j] - Ntip != i + 1 || x2[j] <= Ntip) continue;
		/* can now make the sequence from */
		/* the root to the current node */
		lt = LENGTH(VECTOR_ELT(seqnod, i));
		tmp_vec = allocVector(INTSXP, lt + 1);
		for (k = 0; k < lt; k++)
		  INTEGER(tmp_vec)[k] = INTEGER(VECTOR_ELT(seqnod, i))[k];
		INTEGER(tmp_vec)[lt] = x2[j];
		SET_VECTOR_ELT(seqnod, x2[j] - Ntip - 1, tmp_vec);
	    }
	    done[i] = 1;
	    sumdone++;
	}
    }

    /* build the sequence from root to tip */
    /* by simply looping through 'edge' */
    for (i = 0; i < nbedge; i++) {
        /* skip the internal edges */
        if (x2[i] > Ntip) continue;
	lt = LENGTH(VECTOR_ELT(seqnod, x1[i] - Ntip - 1));
	tmp_vec = allocVector(INTSXP, lt + 1);
	for (j = 0; j < lt; j++)
	  INTEGER(tmp_vec)[j] = INTEGER(VECTOR_ELT(seqnod, x1[i] - Ntip - 1))[j];
	INTEGER(tmp_vec)[lt] = x2[i];
	SET_VECTOR_ELT(ans, x2[i] - 1, tmp_vec);
    }

    UNPROTECT(6);
    return ans;
} /* EOF seq_root2tip */

SEXP bipartition(SEXP edge1, SEXP edge2, SEXP nbtip, SEXP nbnode)
{
    int i, j, k, nbedge, *x1, *x2, lt, lt2, inod, Ntip, Nnode;
    SEXP ans, seqnod, tmp_vec;

    PROTECT(edge1 = coerceVector(edge1, INTSXP));
    PROTECT(edge2 = coerceVector(edge2, INTSXP));
    PROTECT(nbtip = coerceVector(nbtip, INTSXP));
    PROTECT(nbnode = coerceVector(nbnode, INTSXP));
    x1 = INTEGER(edge1);
    x2 = INTEGER(edge2);
    Ntip = *INTEGER(nbtip);
    Nnode = *INTEGER(nbnode);
    nbedge = LENGTH(edge1);

    PROTECT(ans = allocVector(VECSXP, Nnode));

    seqnod = seq_root2tip(edge1, edge2, nbtip, nbnode);

    for (i = 0; i < LENGTH(seqnod); i++) { /* for each tip */
        lt = LENGTH(VECTOR_ELT(seqnod, i));
	for (j = 0; j < lt - 1; j++) {
	    inod = INTEGER(VECTOR_ELT(seqnod, i))[j] - Ntip - 1;
	    if (VECTOR_ELT(ans, inod) == R_NilValue) {
	        tmp_vec = allocVector(INTSXP, 1);
		INTEGER(tmp_vec)[0] = i + 1;
	    } else {
	        lt2 = LENGTH(VECTOR_ELT(ans, inod));
		tmp_vec = allocVector(INTSXP, lt2 + 1);
		for (k = 0; k < lt2; k++)
		  INTEGER(tmp_vec)[k] = INTEGER(VECTOR_ELT(ans, inod))[k];
		INTEGER(tmp_vec)[lt2] = i + 1;
	    }
	    SET_VECTOR_ELT(ans, inod, tmp_vec);
	}
    }

    UNPROTECT(5);
    return ans;
} /* bipartition */
