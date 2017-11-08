/* bipartition.c    2017-07-28 */

/* Copyright 2005-2017 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "ape.h"

SEXP seq_root2tip(SEXP edge, SEXP nbtip, SEXP nbnode)
{
    int i, j, k, Nedge, *x, *done, dn, sumdone, lt, ROOT, Ntip, Nnode;
    SEXP ans, seqnod, tmp_vec;

    /* The following is needed only if we are not sure
       that the storage mode of `edge' is "integer". */
    PROTECT(edge = coerceVector(edge, INTSXP));
    PROTECT(nbtip = coerceVector(nbtip, INTSXP));
    PROTECT(nbnode = coerceVector(nbnode, INTSXP));
    x = INTEGER(edge); /* copy the pointer */
    Ntip = *INTEGER(nbtip);
    Nnode = *INTEGER(nbnode);
    Nedge = LENGTH(edge)/2;
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
	    for (j = 0; j < Nedge; j++) {
	        /* skip the terminal edges, we look only for nodes now */
	        if (x[j] - Ntip != i + 1 || x[j + Nedge] <= Ntip) continue;
		/* can now make the sequence from */
		/* the root to the current node */
		lt = LENGTH(VECTOR_ELT(seqnod, i));
		tmp_vec = allocVector(INTSXP, lt + 1);
		for (k = 0; k < lt; k++)
		  INTEGER(tmp_vec)[k] = INTEGER(VECTOR_ELT(seqnod, i))[k];
		INTEGER(tmp_vec)[lt] = x[j + Nedge];
		SET_VECTOR_ELT(seqnod, x[j + Nedge] - Ntip - 1, tmp_vec);
	    }
	    done[i] = 1;
	    sumdone++;
	}
    }

    /* build the sequence from root to tip */
    /* by simply looping through 'edge' */
    for (i = 0; i < Nedge; i++) {
        /* skip the internal edges */
        if (x[i + Nedge] > Ntip) continue;
	lt = LENGTH(VECTOR_ELT(seqnod, x[i] - Ntip - 1));
	tmp_vec = allocVector(INTSXP, lt + 1);
	for (j = 0; j < lt; j++)
	  INTEGER(tmp_vec)[j] = INTEGER(VECTOR_ELT(seqnod, x[i] - Ntip - 1))[j];
	INTEGER(tmp_vec)[lt] = x[i + Nedge];
	SET_VECTOR_ELT(ans, x[i + Nedge] - 1, tmp_vec);
    }

    UNPROTECT(5);
    return ans;
} /* EOF seq_root2tip */

//SEXP bipartition(SEXP edge, SEXP nbtip, SEXP nbnode)
//{
//    int i, j, k, lt, lt2, inod, Ntip, Nnode;
//    SEXP ans, seqnod, tmp_vec;
//
//    PROTECT(edge = coerceVector(edge, INTSXP));
//    PROTECT(nbtip = coerceVector(nbtip, INTSXP));
//    PROTECT(nbnode = coerceVector(nbnode, INTSXP));
//    Ntip = *INTEGER(nbtip);
//    Nnode = *INTEGER(nbnode);
//
//    PROTECT(ans = allocVector(VECSXP, Nnode));
//    PROTECT(seqnod = seq_root2tip(edge, nbtip, nbnode));
//
//    for (i = 0; i < LENGTH(seqnod); i++) { /* for each tip */
//        lt = LENGTH(VECTOR_ELT(seqnod, i));
//	for (j = 0; j < lt - 1; j++) {
//	    inod = INTEGER(VECTOR_ELT(seqnod, i))[j] - Ntip - 1;
//	    if (VECTOR_ELT(ans, inod) == R_NilValue) {
//	        tmp_vec = allocVector(INTSXP, 1);
//		INTEGER(tmp_vec)[0] = i + 1;
//	    } else {
//	        lt2 = LENGTH(VECTOR_ELT(ans, inod));
//		tmp_vec = allocVector(INTSXP, lt2 + 1);
//		for (k = 0; k < lt2; k++)
//		  INTEGER(tmp_vec)[k] = INTEGER(VECTOR_ELT(ans, inod))[k];
//		INTEGER(tmp_vec)[lt2] = i + 1;
//	    }
//	    SET_VECTOR_ELT(ans, inod, tmp_vec);
//	}
//    }
//
//    UNPROTECT(5);
//    return ans;
//} /* bipartition */

//int SameClade(SEXP clade1, SEXP clade2)
//{
//    int i, n = LENGTH(clade1), *c1, *c2;
//
//    if (n != LENGTH(clade2)) return 0;
//
//    c1 = INTEGER(clade1);
//    c2 = INTEGER(clade2);
//    for (i = 0; i < n; i++)
//      if (c1[i] != c2[i]) return 0;
//
//    return 1;
//}

//SEXP prop_part(SEXP TREES, SEXP nbtree, SEXP keep_partitions)
//{
//    int i, j, k, KeepPartition, Ntree, Ntip, Nnode, Npart, NpartCurrent, *no;
//    SEXP bp, ans, nbtip, nbnode, number;
//
//    PROTECT(nbtree = coerceVector(nbtree, INTSXP));
//    PROTECT(keep_partitions = coerceVector(keep_partitions, INTSXP));
//    Ntree = *INTEGER(nbtree);
//    KeepPartition = *INTEGER(keep_partitions);
//
//
//    Ntip = LENGTH(getListElement(VECTOR_ELT(TREES, 0), "tip.label"));
//    Nnode = *INTEGER(getListElement(VECTOR_ELT(TREES, 0), "Nnode"));
//
//    PROTECT(nbtip = allocVector(INTSXP, 1));
//    PROTECT(nbnode = allocVector(INTSXP, 1));
//    INTEGER(nbtip)[0] = Ntip;
//    INTEGER(nbnode)[0] = Nnode;
//
//    if (KeepPartition) Npart = Ntree * (Ntip - 2) + 1;
//    else Npart = Ntip - 1;
//
//    PROTECT(number = allocVector(INTSXP, Npart));
//    no = INTEGER(number); /* copy the pointer */
//    /* The first partition in the returned list has all tips,
//       so it is observed in all trees: */
//    no[0] = Ntree;
//    /* The partitions in the first tree are obviously observed once: */
//    for (i = 1; i < Nnode; i++) no[i] = 1;
//
//    if (KeepPartition) {
//        for (i = Nnode; i < Npart; i++) no[i] = 0;
//
//        PROTECT(ans = allocVector(VECSXP, Npart));
//	PROTECT(bp = bipartition(getListElement(VECTOR_ELT(TREES, 0), "edge"),
//				 nbtip, nbnode));
//	for (i = 0; i < Nnode; i++)
//	  SET_VECTOR_ELT(ans, i, VECTOR_ELT(bp, i));
//	UNPROTECT(1);
//    } else {
//        PROTECT(ans = bipartition(getListElement(VECTOR_ELT(TREES, 0), "edge"),
//				  nbtip, nbnode));
//    }
//
//    NpartCurrent = Nnode;
//
//    /* We start on the 2nd tree: */
//    for (k = 1; k < Ntree; k++) {
//
///* in case there are trees with multichotomies: */
//	nbnode = getListElement(VECTOR_ELT(TREES, k), "Nnode");
//	Nnode = INTEGER(nbnode)[0];
//
//        PROTECT(bp = bipartition(getListElement(VECTOR_ELT(TREES, k), "edge"),
//				 nbtip, nbnode));
//	for (i = 1; i < Nnode; i++) {
//	    j = 1;
//next_j:
//	    if (SameClade(VECTOR_ELT(bp, i), VECTOR_ELT(ans, j))) {
//	        no[j]++;
//		continue;
//	    }
//	    j++;
//	    if (j < NpartCurrent) goto next_j;
//	    if (KeepPartition) {
//	        no[NpartCurrent]++;
//		SET_VECTOR_ELT(ans, NpartCurrent, VECTOR_ELT(bp, i));
//		NpartCurrent++;
//	    }
//	}
//	UNPROTECT(1);
//    }
//
//    if (KeepPartition && NpartCurrent < Npart) {
//        PROTECT(bp = allocVector(VECSXP, NpartCurrent));
//	for (i = 0; i < NpartCurrent; i++)
//	  SET_VECTOR_ELT(bp, i, VECTOR_ELT(ans, i));
//	setAttrib(bp, install("number"), number);
//	UNPROTECT(7);
//	return bp;
//    } else {
//        setAttrib(ans, install("number"), number);
//	UNPROTECT(6);
//	return ans;
//    }
//} /* prop_part */
