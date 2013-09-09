/* nj.c       2011-10-20 */

/* Copyright 2006-2011 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "ape.h"

double sum_dist_to_i(int n, double *D, int i)
/* returns the sum of all distances D_ij between i and j
   with j = 1...n and j != i */
{
/* we use the fact that the distances are arranged sequentially
   in the lower triangle, e.g. with n = 6 the 15 distances are
   stored as (the C indices are indicated):

           i
     1  2  3  4  5

  2  0
  3  1  5
j 4  2  6  9
  5  3  7 10 12
  6  4  8 11 13 14

  so that we sum the values of the ith column--1st loop--and those of
  (i - 1)th row (labelled 'i')--2nd loop */

	double sum = 0;
	int j, start, end;

	if (i < n) {
		/* the expression below CANNOT be factorized
		   because of the integer operations (it took
		   me a while to find out...) */
		start = n*(i - 1) - i*(i - 1)/2;
		end = start + n - i;
		for (j = start; j < end; j++) sum += D[j];
	}

	if (i > 1) {
		start = i - 2;
		for (j = 1; j <= i - 1; j++) {
			sum += D[start];
			start += n - j - 1;
		}
	}

	return(sum);
}

void C_nj(double *D, int *N, int *edge1, int *edge2, double *edge_length)
{
	double *S, Sdist, Ndist, *new_dist, A, B, smallest_S, x, y;
	int n, i, j, k, ij, smallest, OTU1, OTU2, cur_nod, o_l, *otu_label;

	S = &Sdist;
	new_dist = &Ndist;
	otu_label = &o_l;

	n = *N;
	cur_nod = 2*n - 2;

	S = (double*)R_alloc(n + 1, sizeof(double));
	new_dist = (double*)R_alloc(n*(n - 1)/2, sizeof(double));
	otu_label = (int*)R_alloc(n + 1, sizeof(int));

	for (i = 1; i <= n; i++) otu_label[i] = i; /* otu_label[0] is not used */

	k = 0;

	while (n > 3) {

		for (i = 1; i <= n; i++)
			S[i] = sum_dist_to_i(n, D, i); /* S[0] is not used */

		ij = 0;
		smallest_S = 1e50;
		B = n - 2;
		for (i = 1; i < n; i++) {
			for (j = i + 1; j <= n; j++) {
				A = B*D[ij] - S[i] - S[j];
				if (A < smallest_S) {
					OTU1 = i;
					OTU2 = j;
					smallest_S = A;
					smallest = ij;
				}
				ij++;
			}
		}

		edge2[k] = otu_label[OTU1];
		edge2[k + 1] = otu_label[OTU2];
		edge1[k] = edge1[k + 1] = cur_nod;

		/* get the distances between all OTUs but the 2 selected ones
		   and the latter:
		   a) get the sum for both
		   b) compute the distances for the new OTU */

		A = D[smallest];
		ij = 0;
		for (i = 1; i <= n; i++) {
			if (i == OTU1 || i == OTU2) continue;
			x = D[give_index(i, OTU1, n)]; /* dist between OTU1 and i */
 			y = D[give_index(i, OTU2, n)]; /* dist between OTU2 and i */
			new_dist[ij] = (x + y - A)/2;
			ij++;
		}
		/* compute the branch lengths */
		B = (S[OTU1] - S[OTU2])/B; /* don't need B anymore */
		edge_length[k] = (A + B)/2;
		edge_length[k + 1] = (A - B)/2;

		/* update before the next loop
		   (we are sure that OTU1 < OTU2) */
		if (OTU1 != 1)
			for (i = OTU1; i > 1; i--)
				otu_label[i] = otu_label[i - 1];
		if (OTU2 != n)
			for (i = OTU2; i < n; i++)
				otu_label[i] = otu_label[i + 1];
		otu_label[1] = cur_nod;

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

		cur_nod--;
		k = k + 2;
	}

	for (i = 0; i < 3; i++) {
		edge1[*N*2 - 4 - i] = cur_nod;
		edge2[*N*2 - 4 - i] = otu_label[i + 1];
	}

	edge_length[*N*2 - 4] = (D[0] + D[1] - D[2])/2;
	edge_length[*N*2 - 5] = (D[0] + D[2] - D[1])/2;
	edge_length[*N*2 - 6] = (D[2] + D[1] - D[0])/2;
}
