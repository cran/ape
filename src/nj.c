/* nj.c       2006-01-13 */

/* Copyright 2006 Emmanuel Paradis

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

#define DINDEX(i, j) n*(i - 1) - i*(i - 1)/2 + j - i - 1

int give_index(int i, int j, int n)
{
    if (i > j) return(DINDEX(j, i));
    else return(DINDEX(i, j));
} /* EOF give_index */

double sum_dist_to_i(int n, double *D, int i)
/* returns the sum of all distances D_ij between i and j
   with j between 1, and n and j != i */
{
    double sum;
    int j;

    sum = 0;

    if (i != 1) {
        for (j = 1; j < i; j++)
	  sum += D[DINDEX(j, i)];
    }

    if (i != n) {
        for (j = i + 1; j <= n; j++)
	  sum += D[DINDEX(i, j)];
    }

    return(sum);
} /* EOF sum_dist_to_i */

#define GET_I_AND_J                                               \
/* Find the 'R' indices of the two corresponding OTUs */          \
/* The indices of the first element of the pair in the            \
   distance matrix are n-1 times 1, n-2 times 2, n-3 times 3,     \
   ..., once n-1. Given this, the algorithm below is quite        \
   straightforward.*/                                             \
    i = 0;                                                        \
    for (OTU1 = 1; OTU1 < n; OTU1++) {                            \
        i += n - OTU1;                                            \
	if (i >= smallest + 1) break;                             \
    }                                                             \
    /* Finding the second OTU is easier! */                       \
    OTU2 = smallest + 1 + OTU1 - n*(OTU1 - 1) + OTU1*(OTU1 - 1)/2;

#define SET_CLADE                           \
/* give the node and tip numbers to edge */ \
    edge2[k] = otu_label[OTU1 - 1];         \
    edge2[k + 1] = otu_label[OTU2 - 1];     \
    edge1[k] = edge1[k + 1] = cur_nod;

void nj(double *D, int *N, int *edge1, int *edge2, double *edge_length)
{
    double SUMD, Sdist, *S, Ndist, *new_dist, A, B, *DI, d_i, x, y;
    int n, i, j, k, ij, smallest, OTU1, OTU2, cur_nod, o_l, *otu_label;

    S = &Sdist;
    new_dist = &Ndist;
    otu_label = &o_l;
    DI = &d_i;

    n = *N;
    cur_nod = 2 - n;

    S = (double*)R_alloc(n*(n - 1)/2, sizeof(double));
    new_dist = (double*)R_alloc(n*(n - 1)/2, sizeof(double));
    otu_label = (int*)R_alloc(n, sizeof(int));
    DI = (double*)R_alloc(n - 2, sizeof(double));

    for (i = 0; i < n; i++) otu_label[i] = i + 1;
    k = 0;

    /* First, look if there are any distances equal to 0. */
    /* Since there may be multichotomies, we loop
       through the OTUs instead of the distances. */

    OTU1 = 1;
    while (OTU1 < n) {
        OTU2 = OTU1 + 1;
        while (OTU2 <= n) {
	    if (D[DINDEX(OTU1, OTU2)] == 0) {
		SET_CLADE
		edge_length[k] = edge_length[k + 1] = 0.;
		k = k + 2;

		/* update */

		/* We remove the second tip label: */
		if (OTU2 < n) {
  		    for (i = OTU2; i < n; i++)
		      otu_label[i - 1] = otu_label[i];
		}

		ij = 0;
		for (i = 1; i < n; i++) {
		    if (i == OTU2) continue;
		    for (j = i + 1; j <= n; j++) {
		        if (j == OTU2) continue;
			new_dist[ij] = D[DINDEX(i, j)];
			ij++;
		    }
		}
		n--;
		for (i = 0; i < n*(n - 1)/2; i++) D[i] = new_dist[i];

		otu_label[OTU1 - 1] = cur_nod;
		/* to avoid adjusting the internal branch at the end: */
		DI[-cur_nod - 1] = 0;
		cur_nod++;
	    } else OTU2++;
	}
	OTU1++;
    }

    while (n > 3) {

        SUMD = 0;
	for (i = 0; i < n*(n - 1)/2; i++) SUMD += D[i];

	for (i = 1; i < n; i++) {
	    for (j = i + 1; j <= n; j++) {
	        /* we know that i < j, so: */
	        ij =  DINDEX(i, j);
	        A = sum_dist_to_i(n, D, i) - D[ij];
	        B = sum_dist_to_i(n, D, j) - D[ij];
	        S[ij] = (A + B)/(2*n - 4) + 0.5*D[ij] + (SUMD - A - B - D[ij])/(n - 2);
	    }
	}

	/* find the 'C' index of the smallest value of S */
	smallest = 0;
	for (i = 1; i < n*(n - 1)/2; i++)
	  if (S[smallest] > S[i]) smallest = i;

	GET_I_AND_J

	SET_CLADE

        /* get the distances between all OTUs but the 2 selected ones
           and the latter:
             a) get the sum for both
	     b) compute the distances for the new OTU */
        A = B = ij = 0;
        for (i = 1; i <= n; i++) {
            if (i == OTU1 || i == OTU2) continue;
            x = D[give_index(i, OTU1, n)]; /* distance between OTU1 and i */
            y = D[give_index(i, OTU2, n)]; /* distance between OTU2 and i */
            new_dist[ij] = (x + y)/2;
            A += x;
            B += y;
            ij++;
        }
        /* compute the branch lengths */
        A /= n - 2;
        B /= n - 2;
        edge_length[k] = (D[smallest] + A - B)/2;
        edge_length[k + 1] = (D[smallest] + B - A)/2;
        DI[-cur_nod - 1] = D[smallest];

        /* update before the next loop */
        if (OTU1 > OTU2) { /* make sure that OTU1 < OTU2 */
            i = OTU1;
	    OTU1 = OTU2;
	    OTU2 = i;
        }
        if (OTU1 != 1)
          for (i = OTU1 - 1; i > 0; i--) otu_label[i] = otu_label[i - 1];
        if (OTU2 != n)
          for (i = OTU2; i <= n; i++) otu_label[i - 1] = otu_label[i];
        otu_label[0] = cur_nod;

        for (i = 1; i < n; i++) {
            if (i == OTU1 || i == OTU2) continue;
	    for (j = i + 1; j <= n; j++) {
	        if (j == OTU1 || j == OTU2) continue;
		new_dist[ij] = D[DINDEX(i, j)];
		ij++;
	    }
        }

	n--;
	for (i = 0; i < n*(n - 1)/2; i++) D[i] = new_dist[i];

	cur_nod = cur_nod + 1;
	k = k + 2;
    }

    for (i = 0; i < 3; i++) {
        edge1[*N*2 - 4 - i] = cur_nod;
	edge2[*N*2 - 4 - i] = otu_label[i];
    }

    edge_length[*N*2 - 4] = (D[0] + D[1] - D[2])/2;
    edge_length[*N*2 - 5] = (D[0] + D[2] - D[1])/2;
    edge_length[*N*2 - 6] = (D[2] + D[1] - D[0])/2;

    for (i = 0; i < *N*2 - 3; i++) {
        if (edge2[i] > 0) continue;
	/* In case there are zero branch lengths (see above): */
	if (DI[-edge2[i] - 1] == 0) continue;
	edge_length[i] -= DI[-edge2[i] - 1]/2;
    }
} /* EOF nj */
