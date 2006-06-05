/* bipartition.c    2006-02-01 */

/* Copyright 2005-2006 Emmanuel Paradis */

/* This file is part of the `ape' library for R and related languages. */
/* It is made available under the terms of the GNU General Public */
/* License, version 2, or at your option, any later version, */
/* incorporated herein by reference. */

/* This program is distributed in the hope that it will be */
/* useful, but WITHOUT ANY WARRANTY; without even the implied */
/* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR */
/* PURPOSE.  See the GNU General Public License for more */
/* details. */

/* You should have received a copy of the GNU General Public */
/* License along with this program; if not, write to the Free */
/* Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, */
/* MA 02111-1307, USA */

#include <R.h>
#include <Rinternals.h>

SEXP seq_root2tip(SEXP edge1, SEXP edge2)
{
    int i, j, k, nbedge, nbtip, nbnode, *x1, *x2, *done, dn, sumdone, lt;
    SEXP ans, seqnod, tmp_vec;

    PROTECT(edge1 = coerceVector(edge1, INTSXP));
    PROTECT(edge2 = coerceVector(edge2, INTSXP));
    x1 = INTEGER(edge1);
    x2 = INTEGER(edge2);
    nbedge = LENGTH(edge1);

    /* count the number of tips */
    if (x2[0] > 0) nbtip = 1; else nbtip = 0;
    /* count the number of nodes by finding the smallest element in 'edge1' */
    nbnode = x1[0];
    /* start this loop on the 2nd edge */
    for (i = 1; i < nbedge; i++) {
        if (x2[i] > 0) nbtip++;
	if (x1[i] < nbnode) nbnode = x1[i];
    }
    nbnode = -nbnode;

    PROTECT(ans = allocVector(VECSXP, nbtip));
    PROTECT(seqnod = allocVector(VECSXP, nbnode));

    done = &dn;
    done = (int*)R_alloc(nbnode, sizeof(int));
    for (i = 0; i < nbnode; i++) done[i] = 0;

    tmp_vec = allocVector(INTSXP, 1);
    INTEGER(tmp_vec)[0] = -1;
    SET_VECTOR_ELT(seqnod, 0, tmp_vec);
    sumdone = 0;

    while (sumdone < nbnode) {
        for (i = 0; i < nbnode; i++) { /* loop through all nodes */
	    /* if the vector is not empty and its */
	    /* descendants are not yet found */
	    if (VECTOR_ELT(seqnod, i) == R_NilValue || done[i]) continue;
	    /* look for the descendants in 'edge': */
	    for (j = 0; j < nbedge; j++) {
	        /* skip the terminal edges, we look only for nodes now */
	        if (-x1[j] != i + 1 || x2[j] > 0) continue;
		/* can now make the sequence from */
		/* the root to the current node */
		lt = LENGTH(VECTOR_ELT(seqnod, i));
		tmp_vec = allocVector(INTSXP, lt + 1);
		for (k = 0; k < lt; k++)
		  INTEGER(tmp_vec)[k] = INTEGER(VECTOR_ELT(seqnod, i))[k];
		INTEGER(tmp_vec)[lt] = x2[j];
		SET_VECTOR_ELT(seqnod, -x2[j] - 1, tmp_vec);
	    }
	    done[i] = 1;
	    sumdone++;
	}
    }

    /* build the sequence from root to tip */
    /* by simply looping through 'edge' */
    for (i = 0; i < nbedge; i++) {
        /* skip the internal edges */
        if (x2[i] < 0) continue;
	lt = LENGTH(VECTOR_ELT(seqnod, -x1[i] - 1));
	tmp_vec = allocVector(INTSXP, lt + 1);
	for (j = 0; j < lt; j++)
	  INTEGER(tmp_vec)[j] = INTEGER(VECTOR_ELT(seqnod, -x1[i] - 1))[j];
	INTEGER(tmp_vec)[lt] = x2[i];
	SET_VECTOR_ELT(ans, x2[i] - 1, tmp_vec);
    }

    UNPROTECT(4);
    return(ans);
} /* seq_root2tip */

SEXP bipartition(SEXP edge1, SEXP edge2)
{
    int i, j, k, nbnode, nbedge, *x1, *x2, lt, lt2, inod;
    SEXP ans, seqnod, tmp_vec;

    PROTECT(edge1 = coerceVector(edge1, INTSXP));
    PROTECT(edge2 = coerceVector(edge2, INTSXP));
    x1 = INTEGER(edge1);
    x2 = INTEGER(edge2);
    nbedge = LENGTH(edge1);

    /* count the number of nodes by finding the smallest element in 'edge1' */
    nbnode = x1[0];
    for (i = 1; i < nbedge; i++)
      if (x1[i] < nbnode) nbnode = x1[i];
    nbnode = -nbnode;

    PROTECT(ans = allocVector(VECSXP, nbnode));

    seqnod = seq_root2tip(edge1, edge2);

    for (i = 0; i < LENGTH(seqnod); i++) { /* for each tip */
        lt = LENGTH(VECTOR_ELT(seqnod, i));
	for (j = 0; j < lt - 1; j++) {
	    inod = -1 - INTEGER(VECTOR_ELT(seqnod, i))[j];
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

    UNPROTECT(3);
    return(ans);
} /* bipartition */
