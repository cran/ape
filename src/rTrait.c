/* rTrait.c       2010-05-26 */

/* Copyright 2010 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>

void rTraitCont(int *model, int *Nedge, int *edge1, int *edge2, double *el,
		double *sigma, double *alpha, double *theta, double *x)
{
/* The tree must be in pruningwise order */
	int i;

	switch(*model) {
	case 1 : for (i = *Nedge - 1; i >= 0; i--) {
			GetRNGstate();
			x[edge2[i]] = x[edge1[i]] + sqrt(el[i]) * sigma[i] * norm_rand();
			PutRNGstate();
		}
		break;
	case 2 : for (i = *Nedge - 1; i >= 0; i--) {
			GetRNGstate();
			x[edge2[i]] = x[edge1[i]] + (sigma[i]*norm_rand() - alpha[i]*(x[edge1[i]] - theta[i])) * el[i]; /* need sqrt(el[i]) ? */
			PutRNGstate();
		}
		break;
	}
}
