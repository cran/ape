/* mlphylo.c       2008-01-03 */

/* Copyright 2006-2008 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

typedef struct {
	int *n;
	int *s;
	double *w;
	unsigned char *seq;
	double *anc;
} dna_matrix;

typedef struct {
	int *edge1;
	int *edge2;
	double *el;
} phylo;

typedef struct {
	int *npart;
	int *partition;
	int *model;
	double *xi;
	double *para;
	int *npara;
	double *alpha;
	int *nalpha;
	int *ncat;
	double *invar;
	int *ninvar;
} DNAmodel;

typedef struct {
	dna_matrix X;
	phylo PHY;
	DNAmodel MOD;
	double *BF;
} DNAdata;

typedef struct {
	DNAdata *D; int i;
} info;


void tQ_unbalBF(double *BF, double *P)
/* This function computes the rate matrix Q multiplied by
   time t in the case of unbalanced base frequencies.
   The arguments are:
  BF: the base frequencies
   P: (input) the matrix of substitution rates
      (output) tQ
   NOTE: P must already be multiplied by t */
{
   P[1] *= BF[0];  P[2] *= BF[0];  P[3] *= BF[0];
   P[4] *= BF[1];  P[6] *= BF[1];  P[7] *= BF[1];
   P[8] *= BF[2];  P[9] *= BF[2]; P[11] *= BF[2];
  P[12] *= BF[3]; P[13] *= BF[3]; P[14] *= BF[3];

   P[0] = -P[4] - P[8] - P[12];
   P[5] = -P[1] - P[9] - P[13];
  P[10] = -P[2] - P[6] - P[14];
  P[15] = -P[3] - P[7] - P[11];
}

void mat_expo4x4(double *P)
/* This function computes the exponential of a 4x4 matrix */
{
  double U[16], vl[4], WR[4], Uinv[16], WI[4], work[32];
  int i, info, ipiv[16], n = 4, lw = 32, ord[4];
  char yes = 'V', no = 'N';

  /* The matrix is not symmetric, so we use 'dgeev'. */
  /* We take the real part of the eigenvalues -> WR */
  /* and the right eigenvectors (vr) -> U */
  F77_CALL(dgeev)(&no, &yes, &n, P, &n, WR, WI, vl, &n,
		  U, &n, work, &lw, &info);

  /* It is not necessary to sort the eigenvalues... */
  /* Copy U -> P */
  for (i = 0; i < 16; i++) P[i] = U[i];

  /* For the inversion, we first make Uinv an identity matrix */
  for (i = 1; i < 15; i++) Uinv[i] = 0;
  Uinv[0] = Uinv[5] = Uinv[10] = Uinv[15] = 1;

  /* The matrix is not symmetric, so we use 'dgesv'. */
  /* This subroutine puts the result in Uinv (B) */
  /* (P [= U] is erased) */
  F77_CALL(dgesv)(&n, &n, P, &n, ipiv, Uinv, &n, &info);

  /* The matrix product of U with the eigenvalues diagonal matrix: */
  for (i = 0; i < 4; i++) U[i] *= exp(WR[0]);
  for (i = 4; i < 8; i++) U[i] *= exp(WR[1]);
  for (i = 8; i < 12; i++) U[i] *= exp(WR[2]);
  for (i = 12; i < 16; i++) U[i] *= exp(WR[3]);

  /* The second matrix product with U^-1 */
  P[1] = U[1]*Uinv[0] + U[5]*Uinv[1] + U[9]*Uinv[2] + U[13]*Uinv[3];
  P[2] = U[2]*Uinv[0] + U[6]*Uinv[1] + U[10]*Uinv[2] + U[14]*Uinv[3];
  P[3] = U[3]*Uinv[0] + U[7]*Uinv[1] + U[11]*Uinv[2] + U[15]*Uinv[3];
  P[4] = U[0]*Uinv[4] + U[4]*Uinv[5] + U[8]*Uinv[6] + U[12]*Uinv[7];
  P[6] = U[2]*Uinv[4] + U[6]*Uinv[5] + U[10]*Uinv[6] + U[14]*Uinv[7];
  P[7] = U[3]*Uinv[4] + U[7]*Uinv[5] + U[11]*Uinv[6] + U[15]*Uinv[7];
  P[8] = U[0]*Uinv[8] +  U[4]*Uinv[9] + U[8]*Uinv[10] + U[12]*Uinv[11];
  P[9] = U[1]*Uinv[8] +  U[5]*Uinv[9] + U[9]*Uinv[10] + U[13]*Uinv[11];
  P[11] = U[3]*Uinv[8] +  U[7]*Uinv[9] + U[11]*Uinv[10] + U[15]*Uinv[11];
  P[12] = U[0]*Uinv[12] + U[4]*Uinv[13] + U[8]*Uinv[14] + U[12]*Uinv[15];
  P[13] = U[1]*Uinv[12] + U[5]*Uinv[13] + U[9]*Uinv[14] + U[13]*Uinv[15];
  P[14] = U[2]*Uinv[12] + U[6]*Uinv[13] + U[10]*Uinv[14] + U[14]*Uinv[15];
  P[0] = 1 - P[4] - P[8] - P[12];
  P[5] = 1 - P[1] - P[9] - P[13];
  P[10] = 1 - P[2] - P[6] - P[14];
  P[15] = 1 - P[3] - P[7] - P[11];
}

void PMAT_JC69(double t, double u, double *P)
{
  P[1]=P[2]=P[3]=P[4]=P[6]=P[7]=P[8]=P[9]=P[11]=P[12]=P[13]=P[14]=(1 - exp(-4*u*t))/4;
  P[0] = P[5] = P[10] = P[15] = 1 - 3*P[1];
}

void PMAT_K80(double t, double b, double a, double *P)
{
  double R, p;

  R = a/(2*b);
  p = exp(-2*t/(R + 1));

  P[1] = 0.5*(1 - p); /* A -> C */
  P[2] = 0.25 - 0.5*exp(-t*(2*R + 1)/(R + 1)) + 0.25*p; /* A -> G */
  P[0] = P[5] = P[10] = P[15] = 1 - 2*P[1] - P[2];
  P[3] = P[4] = P[6] = P[11] =  P[9] = P[12] = P[14] = P[1];
  P[7] = P[8] = P[13] = P[2];
}

void PMAT_F81(double t, double u, double *BF, double *P)
{
  double p;
  p = exp(-t*u);

  P[0] = p + (1 - p) * BF[0]; /* A->A */
  P[1] = P[9] = P[13] = (1 - p)*BF[1]; /* A->C, G->C, T->C */
  P[2] = P[6] = P[14] = (1 - p)*BF[2]; /* A->G, C->G, T->G */
  P[3] = P[7] = P[11] = (1 - p)*BF[3]; /* A->T, C->T, G->T */
  P[4] = P[8] = P[12] = (1 - p)*BF[0]; /* C->A, G->A, T->A */
  P[5] = p + P[1]; /* C->C */
  P[10] = p + P[2]; /* G->G */
  P[15] = p + P[3]; /* T->T */
}

void PMAT_F84(double t, double a, double b, double *BF, double *P)
{
  double pI, pII, B, x, y;

  B = exp(-b*t);
  pI = B * (1 - exp(-a*t)); /* probability of at least one event of type I */
  pII = 1 - B; /* probability of at least one event of type II */
  x = pI * (BF[0] + BF[2]);
  y = pI * (BF[1] + BF[3]);

  P[12] = P[4] = pII * BF[0];   /* C->A, T->A */
  P[14] = P[6] = pII * BF[2];   /* C->G, T->G */
  P[9] = P[1] = pII * BF[1];   /* A->C, G->C */
  P[2] = x + P[6];      /* A->G */
  P[11] = P[3] = pII * BF[3];   /* A->T, G->T*/
  P[0] = 1 - P[1] - P[2] - P[3];  /* A->A */
  P[7] = y +  P[3];     /* C->T */
  P[5] = 1 - P[4] - P[6] - P[7];  /* C->C */
  P[8] = x + P[4];      /* G->A */
  P[10] = 1 - P[8] - P[9] - P[11]; /* G->G */
  P[13] = y + P[1];     /* T->C */
  P[15] = 1 - P[12] - P[13] - P[14]; /* T->T */
}

void PMAT_HKY85(double t, double a, double b, double *BF, double *P)
{
  P[2] = P[7] = P[8] = P[13] = t*a;
  P[1] = P[3] = P[4] = P[6] = P[9] = P[11] = P[12] = P[14] = t*b;
  tQ_unbalBF(BF, P);
  mat_expo4x4(P);
}

void PMAT_T92(double t, double a, double b, double *BF, double *P)
{
  double theta, A, B1, B2, C, x, y;

  theta = BF[1] + BF[2];
  A = (1 - theta)/2;
  B1 = (1 + exp(-t));
  B2 = (1 - exp(-t));
  C = exp(-t * ((a/b + 1)/2));
  x = 0.5 * theta * B1;
  y = (1 - theta) * C;

  P[1] = P[6] = P[9] = P[14] = 0.5 * theta * B2; /* A->C, C->G, T->G, G->C */
  P[2] = P[13] = x - theta * C; /* A->G, T->C */
  P[3] = P[4] = P[11] = P[12] = A * B2; /* A->T, C->A, G->T, T->A */
  P[0] = P[15] = 1 - P[1] - P[2] - P[3]; /* A->A, T->T */
  P[5] = P[10] = x + y; /* C->C, G->G */
  P[7] = P[8] = A * B1 - y; /* C->T, G->A */
}

void PMAT_TN93(double t, double a1, double a2, double b,
	       double *BF, double *P)
{

  double A1, A2, B;

  A1 = (1 - exp(-a1*t)); /* transitions between purines (A <-> G) */
  A2 = (1 - exp(-a2*t)); /* transitions between pyrimidines (C <-> T) */
  B = (1 - exp(-b*t));

  P[1] = B * BF[1];                  /* A->C */
  P[2] = A1 * BF[2];                 /* A->G */
  P[3] = B * BF[3];                  /* A->T */
  P[0] = 1 - P[1] - P[2] - P[3];     /* A->A */
  P[4] = B * BF[0];                  /* C->A */
  P[6] = B * BF[2];                  /* C->G */
  P[7] = A2 * BF[3] ;                /* C->T */
  P[5] = 1 - P[4] - P[6] - P[7];     /* C->C */
  P[8] = A1 * BF[0];                 /* G->A */
  P[9] = B * BF[1];                  /* G->C */
  P[11] = B * BF[3];                 /* G->T */
  P[10] = 1 - P[8] - P[9] - P[11];   /* G->G */
  P[12] = B * BF[0];                 /* T->A */
  P[13] = A2 * BF[1];                /* T->C */
  P[14] = B * BF[2];                 /* T->G */
  P[15] = 1 - P[12] - P[13] - P[14]; /* T->T */
}

void PMAT_GTR(double t, double a, double b, double c, double d, double e,
	      double f, double *BF, double *P)
{
  P[1] = P[4] = t*a;
  P[2] = P[8] = t*b;
  P[3] = P[12] = t*c;
  P[6] = P[9] = t*d;
  P[7] = P[13] = t*e;
  P[11] = P[14] = t*f;
  tQ_unbalBF(BF, P);
  mat_expo4x4(P);
}

#define GET_DNA_PARAMETERS \
    /* get the substitution parameters */ \
    model = *(D->MOD.model); \
    /* If the model is not JC69 or F81: */ \
    if (model != 1 && model != 3) { \
        for(i = 0; i < *(D->MOD.npara); i++) \
            u[i] = D->MOD.para[i]; \
    } \
    /* get the shape parameter and calculate the coefficients */ \
    ncat = *(D->MOD.ncat); \
    if (ncat > 1) { \
        /* use `tmp[0]' to store the mean of the coefficients */ \
        /* in order to rescale them */ \
        tmp[0] = 0.; \
        if (*(D->MOD.nalpha) > 1) alpha = *(D->MOD.alpha); \
        else alpha = D->MOD.alpha[k]; \
	for (j = 0; j < ncat; j++) { \
		coef_gamma[j] = qgamma((0.5 + j)/ncat, alpha, \
				       1/alpha, 1, 0); \
                tmp[0] += coef_gamma[j]; \
	} \
        tmp[0] /= ncat; \
        for (j = 0; j < ncat; j++) \
          coef_gamma[j] /= tmp[0]; \
    } else coef_gamma[0] = 1.; \
    /* get the proportion of invariants */ \
    if (*(D->MOD.ninvar)) { \
        if (*(D->MOD.ninvar) > 1) I = *(D->MOD.invar); \
        else I = D->MOD.invar[k]; \
    } else I = 0.; \

void getSiteLik(int n, int d, int j, int nr, DNAdata *D, double *L)
{
	int i;

	if (d <= n - 1) {
		i = d + j*n;
		memset(L, 0, 4*sizeof(double));
		if (D->X.seq[i] & 128) L[0] = 1;
		if (D->X.seq[i] & 64) L[1] = 1;
		if (D->X.seq[i] & 32) L[2] = 1;
		if (D->X.seq[i] & 16) L[3] = 1;
	} else {
		i = (d - n) + j*(n - 2);
		L[0] = D->X.anc[i];
		L[1] = D->X.anc[i + nr];
		L[2] = D->X.anc[i + 2*nr];
		L[3] = D->X.anc[i + 3*nr];
	}
}

#define LOOP_THROUGH_SITES \
    for(j = start; j < end; j++) { \
        memset(tmp, 0, 4*sizeof(double)); \
        getSiteLik(n, d1, j, nr, D, L1); \
        getSiteLik(n, d2, j, nr, D, L2); \
	for(i = 0; i < ncat; i++) { \
	    switch(model) { \
	    case 1 : PMAT_JC69(l1, coef_gamma[i], P1); \
                     PMAT_JC69(l2, coef_gamma[i], P2); break; \
	    case 2 : PMAT_K80(l1, coef_gamma[i], u[0], P1); \
                     PMAT_K80(l2, coef_gamma[i], u[0], P2); break; \
	    case 3 : PMAT_F81(l1, coef_gamma[i], BF, P1); \
                     PMAT_F81(l2, coef_gamma[i], BF, P2); break; \
	    case 4 : PMAT_F84(l1, coef_gamma[i], u[0], BF, P1); \
                     PMAT_F84(l2, coef_gamma[i], u[0], BF, P2); break; \
	    case 5 : PMAT_HKY85(l1, coef_gamma[i], u[0], BF, P1); \
                     PMAT_HKY85(l2, coef_gamma[i], u[0], BF, P2); break; \
	    case 6 : PMAT_T92(l1, coef_gamma[i], u[0], BF, P1); \
                     PMAT_T92(l2, coef_gamma[i], u[0], BF, P2); break; \
	    case 7 : PMAT_TN93(l1, coef_gamma[i], u[0], u[1], BF, P1); \
                     PMAT_TN93(l2, coef_gamma[i], u[0], u[1], BF, P2); break; \
	    case 8 : PMAT_GTR(l1, coef_gamma[i], u[0], u[1], u[2], u[3], u[4], BF, P1); \
                     PMAT_GTR(l2, coef_gamma[i], u[0], u[1], u[2], u[3], u[4], BF, P2); break; \
	    } \
            tmp[0] += (L1[0]*P1[0] + L1[1]*P1[1] + L1[2]*P1[2] + L1[3]*P1[3]) * \
		      (L2[0]*P2[0] + L2[1]*P2[1] + L2[2]*P2[2] + L2[3]*P2[3]); \
            tmp[1] += (L1[0]*P1[4] + L1[1]*P1[5] + L1[2]*P1[6] + L1[3]*P1[7]) * \
		      (L2[0]*P2[4] + L2[1]*P2[5] + L2[2]*P2[6] + L2[3]*P2[7]); \
            tmp[2] += (L1[0]*P1[8] + L1[1]*P1[9] + L1[2]*P1[10] + L1[3]*P1[11]) * \
		      (L2[0]*P2[8] + L2[1]*P2[9] + L2[2]*P2[10] + L2[3]*P2[11]); \
            tmp[3] += (L1[0]*P1[12] + L1[1]*P1[13] + L1[2]*P1[14] + L1[3]*P1[15]) * \
		      (L2[0]*P2[12] + L2[1]*P2[13] + L2[2]*P2[14] + L2[3]*P2[15]); \
	} \
        if (ncat > 1) { \
            tmp[0] /= ncat; \
            tmp[1] /= ncat; \
            tmp[2] /= ncat; \
            tmp[3] /= ncat; \
        } \
	if (D->MOD.ninvar) { \
	    V = 1. - I; \
            tmp[0] = V*tmp[0] + I*L1[0]*L2[0]; \
            tmp[1] = V*tmp[1] + I*L1[1]*L2[1]; \
            tmp[2] = V*tmp[2] + I*L1[2]*L2[2]; \
            tmp[3] = V*tmp[3] + I*L1[3]*L2[3]; \
	} \
        ind = anc - n + j*(n - 2); \
        D->X.anc[ind] = tmp[0]; \
        D->X.anc[ind + nr] = tmp[1]; \
        D->X.anc[ind + 2*nr] = tmp[2]; \
        D->X.anc[ind + 3*nr] = tmp[3]; \
    }

void lik_dna_node(DNAdata *D, int ie)
/*
This function computes the likelihoods at a node for all
nucleotides.
*/
{
	int d1, d2, anc, n, nr;
	int i, j, k, start, end, ind, i1, i2, ncat, model;
	double tmp[4], L1[4], L2[4], l1, l2, P1[16], P2[16], V, coef_gamma[10], u[6], I, alpha, *BF;

	n = *(D->X.n);
	nr = *(D->X.s) * (n - 2);
	BF = D->BF;

	d1 = D->PHY.edge2[ie];
	d2 = D->PHY.edge2[ie + 1];
	anc = D->PHY.edge1[ie];

	/* change these to use them as indices */
	d1--; d2--; anc--;

	l1 = D->PHY.el[ie];
	l2 = D->PHY.el[ie + 1];

	for(k = 0; k < *(D->MOD.npart); k++) {
		start = D->MOD.partition[k*2] - 1;
		end = D->MOD.partition[k*2 + 1] - 1;

		GET_DNA_PARAMETERS

		if (k > 0) {
			l1 *= D->MOD.xi[k - 1];
			l2 *= D->MOD.xi[k - 1];
		}

		LOOP_THROUGH_SITES
	}
} /* lik_dna_node */

void lik_dna_root(DNAdata *D)
/*
This function computes the likelihoods at the root for all
nucleotides.
*/
{
	int i, j, k, start, end, ind, ncat, model, d1, d2, d3, n, N, nr;
	double tmp[4],  L1[4], L2[4], L3[4], l1, l2, l3, P1[16], P2[16], P3[16], V, coef_gamma[10], u[6], I, alpha, *BF;

	n = *(D->X.n); /* number of tips */
	N = 2*n - 3; /* number of edges */
	nr = *(D->X.s) * (n - 2);
	BF = D->BF;

	d1 = D->PHY.edge2[N - 3];
	d2 = D->PHY.edge2[N - 2];
	d3 = D->PHY.edge2[N - 1];

	/* change these to use them as indices */
	d1--; d2--; d3--;

	l1 = D->PHY.el[N - 3];
	l2 = D->PHY.el[N - 2];
	l3 = D->PHY.el[N - 1];

	for(k = 0; k < *(D->MOD.npart); k++) {
		start = D->MOD.partition[k*2] - 1;
		end = D->MOD.partition[k*2 + 1] - 1;

		GET_DNA_PARAMETERS

		if (k > 0) {
			l1 *= D->MOD.xi[k - 1];
			l2 *= D->MOD.xi[k - 1];
			l3 *= D->MOD.xi[k - 1];
		}

		for(j = start; j < end; j++) {
			getSiteLik(n, d1, j, nr, D, L1);
			getSiteLik(n, d2, j, nr, D, L2);
			getSiteLik(n, d3, j, nr, D, L3);
			memset(tmp, 0, 4*sizeof(double));
			for(i = 0; i < ncat; i++) {
				switch(model) {
				case 1 : PMAT_JC69(l1, coef_gamma[i], P1);
					PMAT_JC69(l2, coef_gamma[i], P2);
					PMAT_JC69(l3, coef_gamma[i], P3); break;
				case 2 : PMAT_K80(l1, coef_gamma[i], u[0], P1);
					PMAT_K80(l2, coef_gamma[i], u[0], P2);
					PMAT_K80(l3, coef_gamma[i], u[0], P3); break;
				case 3 : PMAT_F81(l1, coef_gamma[i], BF, P1);
					PMAT_F81(l2, coef_gamma[i], BF, P3);
					PMAT_F81(l3, coef_gamma[i], BF, P3); break;
				case 4 : PMAT_F84(l1, coef_gamma[i], u[0], BF, P1);
					PMAT_F84(l2, coef_gamma[i], u[0], BF, P2);
					PMAT_F84(l3, coef_gamma[i], u[0], BF, P3); break;
				case 5 : PMAT_HKY85(l1, coef_gamma[i], u[0], BF, P1);
					PMAT_HKY85(l2, coef_gamma[i], u[0], BF, P2);
					PMAT_HKY85(l3, coef_gamma[i], u[0], BF, P3); break;
				case 6 : PMAT_T92(l1, coef_gamma[i], u[0], BF, P1);
					PMAT_T92(l2, coef_gamma[i], u[0], BF, P2);
					PMAT_T92(l3, coef_gamma[i], u[0], BF, P3); break;
				case 7 : PMAT_TN93(l1, coef_gamma[i], u[0], u[1], BF, P1);
					PMAT_TN93(l2, coef_gamma[i], u[0], u[1], BF, P2);
					PMAT_TN93(l3, coef_gamma[i], u[0], u[1], BF, P3);break;
				case 8 : PMAT_GTR(l1, coef_gamma[i], u[0], u[1], u[2], u[3], u[4], BF, P1);
					PMAT_GTR(l2, coef_gamma[i], u[0], u[1], u[2], u[3], u[4], BF, P2);
					PMAT_GTR(l3, coef_gamma[i], u[0], u[1], u[2], u[3], u[4], BF, P3); break;
				}
				tmp[0] += (L1[0]*P1[0] + L1[1]*P1[1] + L1[2]*P1[2] + L1[3]*P1[3]) *
					(L2[0]*P2[0] + L2[1]*P2[1] + L2[2]*P2[2] + L2[3]*P2[3]) *
					(L3[0]*P3[0] + L3[1]*P3[1] + L3[2]*P3[2] + L3[3]*P3[3]);
				tmp[1] += (L1[0]*P1[4] + L1[1]*P1[5] + L1[2]*P1[6] + L1[3]*P1[7]) *
					(L2[0]*P2[4] + L2[1]*P2[5] + L2[2]*P2[6] + L2[3]*P2[7]) *
					(L3[0]*P3[4] + L3[1]*P3[5] + L3[2]*P3[6] + L3[3]*P3[7]);
				tmp[2] += (L1[0]*P1[8] + L1[1]*P1[9] + L1[2]*P1[10] + L1[3]*P1[11]) *
					(L2[0]*P2[8] + L2[1]*P2[9] + L2[2]*P2[10] + L2[3]*P2[11]) *
					(L3[0]*P3[8] + L3[1]*P3[9] + L3[2]*P3[10] + L3[3]*P3[11]);
				tmp[3] += (L1[0]*P1[12] + L1[1]*P1[13] + L1[2]*P1[14] + L1[3]*P1[15]) *
					(L2[0]*P2[12] + L2[1]*P2[13] + L2[2]*P2[14] + L2[3]*P2[15]) *
					(L3[0]*P3[12] + L3[1]*P3[13] + L3[2]*P3[14] + L3[3]*P3[15]);
			}
			if (ncat > 1) {
				tmp[0] /= ncat;
				tmp[1] /= ncat;
				tmp[2] /= ncat;
				tmp[3] /= ncat;
			}
			if (D->MOD.ninvar) {
				V = 1. - I;
				tmp[0] = V*tmp[0] + I*L1[0]*L2[0]*L3[0];
				tmp[1] = V*tmp[1] + I*L1[1]*L2[1]*L3[1];
				tmp[2] = V*tmp[2] + I*L1[2]*L2[2]*L3[2];
				tmp[3] = V*tmp[3] + I*L1[3]*L2[3]*L3[3];
			}
			ind = j*(n - 2);
			D->X.anc[ind] = tmp[0];
			D->X.anc[ind + nr] = tmp[1];
			D->X.anc[ind + 2*nr] = tmp[2];
			D->X.anc[ind + 3*nr] = tmp[3];
		}
	}
} /* lik_dna_root */

void lik_dna_tree(DNAdata *D, double *loglik)
{
    int i, j, n, nnode, nsite, nr;
    double tmp;

    n = *(D->X.n);
    nnode = n - 2;
    nsite = *(D->X.s);
    nr = nsite*nnode;

    /* initialize before looping through the tree */
    memset(D->X.anc, 1., nr*4*sizeof(double));

    /* loop through the tree
       We don't do the root node here, so i goes between 0 and 2n - 6 */
    for(i = 0; i < 2*n - 6; i += 2)
	    lik_dna_node(D, i);

    /* We now do the root */
    lik_dna_root(D);
    *loglik = 0.;
    for(j = 0; j < nsite; j++) {
	    tmp = 0.;
	    for (i = 0; i < 4; i++)
		    tmp += D->BF[i] * D->X.anc[j + i*nr];
	    *loglik += D->X.w[j]*log(tmp);
    }
} /* lik_dna_tree */

double fcn_mlphylo_invar(double I, info *INFO)
{
    double loglik;

    INFO->D->MOD.invar[INFO->i] = I;
    lik_dna_tree(INFO->D, &loglik);

    return -loglik;
}

void mlphylo_invar(int N, DNAdata *D, double *loglik)
/*
optimize proportion of invariants
*/
{
    int i;
    info INFO, *infptr;
    double I;

    infptr = &INFO;
    INFO.D = D;

    for(i = 0; i < N; i++) {
        infptr->i = i;
	I = Brent_fmin(0.0, 1.0,
		       (double (*)(double, void*)) fcn_mlphylo_invar,
		       infptr, 1.e-9);
	D->MOD.invar[i] = I;
    }
}

double fcn_mlphylo_gamma(double a, info *INFO)
{
    double loglik;

    INFO->D->MOD.alpha[INFO->i] = a;
    lik_dna_tree(INFO->D, &loglik);

    return -loglik;
}

void mlphylo_gamma(int N, DNAdata *D, double *loglik)
/*
optimize gamma (ISV) parameters
*/
{
    int i;
    info INFO, *infptr;
    double a;

    infptr = &INFO;
    INFO.D = D;

    for(i = 0; i < N; i++) {
        infptr->i = i;
	a = Brent_fmin(0.0, 1.e4,
		       (double (*)(double, void*)) fcn_mlphylo_gamma,
		       infptr, 1.e-6);
	D->MOD.alpha[i] = a;
    }
}

double fcn_mlphylo_para(double p, info *INFO)
{
    double loglik;

    INFO->D->MOD.para[INFO->i] = p;
    lik_dna_tree(INFO->D, &loglik);

    return -loglik;
}

void mlphylo_para(int N, DNAdata *D, double *loglik)
/*
optimize the contrast parameter(s) xi
*/
{
    int i;
    info INFO, *infptr;
    double p;

    infptr = &INFO;
    INFO.D = D;

    for(i = 0; i < N; i++) {
        infptr->i = i;
	p = Brent_fmin(0, 1.e3,
		       (double (*)(double, void*)) fcn_mlphylo_para,
		       infptr, 1.e-6);
	D->MOD.para[i] = p;
    }
}

double fcn_mlphylo_xi(double XI, info *INFO)
{
    double loglik;

    INFO->D->MOD.xi[INFO->i] = XI;
    lik_dna_tree(INFO->D, &loglik);

    return -loglik;
}

void mlphylo_xi(int N, DNAdata *D, double *loglik)
/*
optimize the contrast parameter(s) xi
*/
{
    int i;
    info INFO, *infptr;
    double XI;

    infptr = &INFO;
    INFO.D = D;

    /* In the following, the range of the search algo was changed from */
    /* 0-1000 to 0-100 to avoid infinite looping. (2006-04-15) */
    /* This was changed again to 0-20. (2006-07-17) */

    for(i = 0; i < N; i++) {
        infptr->i = i;
	XI = Brent_fmin(0.0, 2.e1,
		       (double (*)(double, void*)) fcn_mlphylo_xi,
		       infptr, 1.e-4);
	D->MOD.xi[i] = XI;
    }
}

double fcn_mlphylo_edgelength(double l, info *INFO)
{
    double loglik;

    INFO->D->PHY.el[INFO->i] = l;
    lik_dna_tree(INFO->D, &loglik);

    return -loglik;
}

void mlphylo_edgelength(int N, DNAdata *D, double *loglik)
/*
optimize branch lengths
*/
{
    int i;
    info INFO, *infptr;
    double l;

    infptr = &INFO;
    INFO.D = D;

    for(i = 0; i < N; i++) {
        infptr->i = i;
	l = Brent_fmin(0.0, 0.1,
		       (double (*)(double, void*)) fcn_mlphylo_edgelength,
		       infptr, 1.e-6);
	D->PHY.el[i] = l;
    }
}

void mlphylo_DNAmodel(int *n, int *s, unsigned char *SEQ, double *ANC,
		      double *w, int *edge1, int *edge2,
		      double *edge_length, int *npart, int *partition,
		      int *model, double *xi, double *para, int *npara,
		      double *alpha, int *nalpha, int *ncat,
		      double *invar, int *ninvar, double *BF,
		      int *search_tree, int *fixed, double *loglik)
/*
This function iterates to find the MLEs of the substitution
paramaters and of the branch lengths for a given tree.
*/
{
	DNAdata *D, data;

	D = &data;

	D->X.n = n;
	D->X.s = s;
	D->X.w = w;
	D->X.seq = SEQ;
	D->X.anc = ANC;

	D->PHY.edge1 = edge1;
	D->PHY.edge2 = edge2;
	D->PHY.el = edge_length;

	D->MOD.npart = npart;
	D->MOD.partition = partition;
	D->MOD.model = model;
	D->MOD.xi = xi;
	D->MOD.para = para;
	D->MOD.npara = npara;
	D->MOD.alpha = alpha;
	D->MOD.nalpha = nalpha;
	D->MOD.ncat = ncat;
	D->MOD.invar = invar;
	D->MOD.ninvar = ninvar;

	D->BF = BF;

	lik_dna_tree(D, loglik);
	if (! *fixed) {
		if (*npart > 1) mlphylo_xi(*npart - 1, D, loglik);
		if (*npara) mlphylo_para(*npara, D, loglik);
		if (*nalpha) mlphylo_gamma(*nalpha, D, loglik);
		if (*ninvar) mlphylo_invar(*ninvar, D, loglik);
	}
	lik_dna_tree(D, loglik);
} /* mlphylo_DNAmodel */
/*
void jc69(double *P, double *t, double *u)
{
	PMAT_JC69(*t, *u, P);
}

void k80(double *P, double *t, double *u)
{
	PMAT_K80(*t, u[0], u[1], P);
}
*/
