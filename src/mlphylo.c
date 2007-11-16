/* mlphylo.c       2007-09-27 */

/* Copyright 2006-2007 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

typedef struct {
	double *A, *C, *G, *T;
} dna_matrix;

typedef struct {
	int *s1, *s2, *a;
} matching;

typedef struct {
	int *n; int *limit; int *model; double *xi;
} part_dna;

typedef struct {
	double *para; int *n; int *pim;
} para_dna;

typedef struct {
	int *ncat; double *alpha; int *n; int *pim;
} gamma_dna;

typedef struct {
	double *I; int *n; int *pim;
} invar_dna;

typedef struct {
	int *n; int *s; dna_matrix X; double *w;
	matching match; double *edge_length;
	part_dna partition; para_dna PAR;
	gamma_dna GAMMA; invar_dna INV; double *BF;
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

void PMAT_K80(double t, double a, double b, double *P)
{
  double R, p;

  R = a/(2*b);
  p = exp(-2*t/(R + 1));

  P[1] = 0.5 * (1 - p); /* A -> C */
  P[2] = 0.25 - 0.5 * exp(-t*(2*R + 1)/(R + 1)) + 0.25 * p; /* A -> G */
  P[0] = P[5] = P[10] = P[15] = 1 - 2 * P[1] - P[2];
  P[3] = P[4] = P[6] = P[8] =  P[9] = P[12] = P[13] = P[1];
  P[7] = P[11] = P[14] = P[2];
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

#define GET_DNA_PARAMETERS                                        \
    /* get the substitution parameters */                         \
    model = PART.model[k];                                        \
    /* If the model is not JC69 or F81: */                        \
    if (model != 1 && model != 3) {                               \
        ind = 0;                                                  \
        for(i = 0; i < *(PARA.n); i++) {                          \
            if (PARA.pim[i + *(PARA.n)*k]) {                      \
        	    u[ind] = PARA.para[i];                        \
        	    ind++;                                        \
	    }                                                     \
        }                                                         \
    }                                                             \
    /* get the shape parameter and calculate the coefficients */  \
    ncat = 1;                                                     \
    if (GAMMA.ncat[k] > 1) {                                      \
        ncat = GAMMA.ncat[k];                                     \
        /* use `tmp[0]' to store the mean of the coefficients */  \
        /* in order to rescale them */                            \
        tmp[0] = 0.;                                              \
	for(i = 0; i < *(GAMMA.n); i++) {                         \
	    if (GAMMA.pim[i + *(GAMMA.n)*k]) {                    \
	        for (j = 0; j < ncat; j++) {                      \
		    coef_gamma[j] = qgamma((0.5 + j)/ncat,        \
                                           GAMMA.alpha[i],        \
					   1/GAMMA.alpha[i],      \
					   1, 0);                 \
                    tmp[0] += coef_gamma[j];                      \
		}                                                 \
                tmp[0] /= ncat;                                   \
                for (j = 0; j < ncat; j++)                        \
                  coef_gamma[j] /= tmp[0];                        \
	    }                                                     \
	}                                                         \
    } else coef_gamma[0] = 1.;                                    \
    /* get the proportion of invariants */                        \
    I = 0;                                                        \
    for(i = 0; i < *(INV.n); i++) {                               \
        if (INV.pim[i + *(INV.n)*k]) {                            \
	    I = INV.I[i];                                         \
	    break;                                                \
	}                                                         \
    }

#define LOOP_THROUGH_SITES \
    for(j = start; j < end; j++) { \
        i1 = d1 + j*nr; \
        i2 = d2 + j*nr; \
        tmp[0] = tmp[1] = tmp[2] = tmp[3] = 0.0; \
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
	    case 8 : PMAT_GTR(l1, coef_gamma[i], u[0], u[1], \
			      u[2], u[3], u[4], BF, P1); \
PMAT_GTR(l2, coef_gamma[i], u[0], u[1], \
			      u[2], u[3], u[4], BF, P2); break; \
	    } \
	    tmp[0] += (X.A[i1] * P1[0] + X.C[i1] * P1[1] + \
		      X.G[i1] * P1[2] + X.T[i1] * P1[3]) * \
                     (X.A[i2] * P2[0] + X.C[i2] * P2[1] + \
		      X.G[i2] * P2[2] + X.T[i2] * P2[3]); \
	    tmp[1] += (X.A[i1] * P1[4] + X.C[i1] * P1[5] + \
		      X.G[i1] * P1[6] + X.T[i1] * P1[7]) * \
                     (X.A[i2] * P2[4] + X.C[i2] * P2[5] + \
		      X.G[i2] * P2[6] + X.T[i2] * P2[7]); \
	    tmp[2] += (X.A[i1] * P1[8] + X.C[i1] * P1[9] + \
		      X.G[i1] * P1[10] + X.T[i1] * P1[11]) * \
                     (X.A[i2] * P2[8] + X.C[i2] * P2[9] + \
		      X.G[i2] * P2[10] + X.T[i2] * P2[11]); \
	    tmp[3] += (X.A[i1] * P1[12] + X.C[i1] * P1[13] + \
		      X.G[i1] * P1[14] + X.T[i1] * P1[15]) * \
                     (X.A[i2] * P2[12] + X.C[i2] * P2[13] + \
		      X.G[i2] * P2[14] + X.T[i2] * P2[15]); \
	} \
        if (ncat > 1) { \
            tmp[0] /= ncat; \
            tmp[1] /= ncat; \
            tmp[2] /= ncat; \
            tmp[3] /= ncat; \
        } \
        ind = anc + j*nr; \
	if (I > 0) { \
	    V = 1. - I; \
            tmp[0] = V*tmp[0] + I*X.A[i1]*X.A[i2]; \
            tmp[1] = V*tmp[1] + I*X.C[i1]*X.C[i2]; \
            tmp[2] = V*tmp[2] + I*X.G[i1]*X.G[i2]; \
            tmp[3] = V*tmp[3] + I*X.T[i1]*X.T[i2]; \
	} \
        X.A[ind] = tmp[0]; \
        X.C[ind] = tmp[1]; \
        X.G[ind] = tmp[2]; \
        X.T[ind] = tmp[3]; \
Rprintf("");\
    }

/* <FIXME>
The Rprintf() above is needed to avoid a memory corruption,
but I don't know why! (2007-03-27)
</FIXME> */

void lik_dna_node(dna_matrix X, int d1, int d2, int anc,
		  double *edge_length, part_dna PART,
		  para_dna PARA, gamma_dna GAMMA,
		  invar_dna INV, double *BF, int nr)
/*
This function computes the likelihoods at a node for all
nucleotides.
*/
{
    int i, j, k, start, end, ind, i1, i2, ncat, model;
    double tmp[4], l1, l2, P1[16], P2[16], V, coef_gamma[10], u[6], I;

    /* change d1, d2, and anc to use them as indices */
    d1--; d2--; anc--;

    l1 = edge_length[d1];
    l2 = edge_length[d2];

    for(k = 0; k < *(PART.n); k++) {
        start = PART.limit[k*2] - 1;
        end = PART.limit[k*2 + 1];

        GET_DNA_PARAMETERS

	if (k > 0) {
	    l1 *= PART.xi[k - 1];
	    l2 *= PART.xi[k - 1];
	}

	LOOP_THROUGH_SITES
    }
} /* EOF lik_dna_node */

void lik_dna_root(dna_matrix X, int d1, int d2, int d3, int anc,
		  double *edge_length, part_dna PART,
		  para_dna PARA, gamma_dna GAMMA,
		  invar_dna INV, double *BF, int nr)
/*
This function computes the likelihoods at the root for all
nucleotides.
*/
{
    int i, j, k, start, end, ind, ncat, model, i1, i2, i3;
    double tmp[4], l1, l2, l3, P1[16], P2[16], P3[16], V,
      coef_gamma[10], u[6], I;

    /* change d1, d2, d3, and anc to use them as indices */
    d1--; d2--; d3--; anc--;

    l1 = edge_length[d1];
    l2 = edge_length[d2];
    l3 = edge_length[d3];

    for(k = 0; k < *(PART.n); k++) {
        start = PART.limit[k*2] - 1;
        end = PART.limit[k*2 + 1];

        GET_DNA_PARAMETERS

	if (k > 0) {
	    l1 *= PART.xi[k - 1];
	    l2 *= PART.xi[k - 1];
	    l3 *= PART.xi[k - 1];
	}

	for(j = start; j < end; j++) {
	    i1 = d1 + j*nr;
	    i2 = d2 + j*nr;
	    i3 = d3 + j*nr;
	    tmp[0] = tmp[1] = tmp[2] = tmp[3] = 0.;
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
		case 8 : PMAT_GTR(l1, coef_gamma[i], u[0], u[1],
				  u[2], u[3], u[4], BF, P1);
		  PMAT_GTR(l2, coef_gamma[i], u[0], u[1],
			   u[2], u[3], u[4], BF, P2);
		  PMAT_GTR(l3, coef_gamma[i], u[0], u[1],
			   u[2], u[3], u[4], BF, P3); break;
		}
		tmp[0] += (X.A[i1] * P1[0] + X.C[i1] * P1[1] +
			   X.G[i1] * P1[2] + X.T[i1] * P1[3]) *
		  (X.A[i2] * P2[0] + X.C[i2] * P2[1] +
		   X.G[i2] * P2[2] + X.T[i2] * P2[3]) *
		  (X.A[i3] * P3[0] + X.C[i3] * P3[1] +
		   X.G[i3] * P3[2] + X.T[i3] * P3[3]);
		tmp[1] += (X.A[i1] * P1[4] + X.C[i1] * P1[5] +
			   X.G[i1] * P1[6] + X.T[i1] * P1[7]) *
		  (X.A[i2] * P2[4] + X.C[i2] * P2[5] +
		   X.G[i2] * P2[6] + X.T[i2] * P2[7]) *
		  (X.A[i3] * P3[4] + X.C[i3] * P3[5] +
		   X.G[i3] * P3[6] + X.T[i3] * P3[7]);
		tmp[2] += (X.A[i1] * P1[8] + X.C[i1] * P1[9] +
			   X.G[i1] * P1[10] + X.T[i1] * P1[11]) *
		  (X.A[i2] * P2[8] + X.C[i2] * P2[9] +
		   X.G[i2] * P2[10] + X.T[i2] * P2[11]) *
		  (X.A[i3] * P3[8] + X.C[i3] * P3[9] +
		   X.G[i3] * P3[10] + X.T[i3] * P3[11]);
		tmp[3] += (X.A[i1] * P1[12] + X.C[i1] * P1[13] +
			   X.G[i1] * P1[14] + X.T[i1] * P1[15]) *
		  (X.A[i2] * P2[12] + X.C[i2] * P2[13] +
		   X.G[i2] * P2[14] + X.T[i2] * P2[15]) *
		  (X.A[i3] * P3[12] + X.C[i3] * P3[13] +
		   X.G[i3] * P3[14] + X.T[i3] * P3[15]);
	    }
	    if (ncat > 1) {
	        tmp[0] /= ncat;
		tmp[1] /= ncat;
		tmp[2] /= ncat;
		tmp[3] /= ncat;
	    }
	    if (I > 0) { /* maybe use a better comparison */
	        V = 1. - I;
		tmp[0] = V * tmp[0] + I * X.A[ind];
		tmp[1] = V * tmp[1] + I * X.C[ind];
		tmp[2] = V * tmp[2] + I * X.G[ind];
		tmp[3] = V * tmp[3] + I * X.T[ind];
	    }
	    ind = anc + j*nr;
	    X.A[ind] = tmp[0];
	    X.C[ind] = tmp[1];
	    X.G[ind] = tmp[2];
	    X.T[ind] = tmp[3];
	}
    }
} /* EOF lik_dna_root */

/*----------------------------------------------------*\
| Il faudra ajouter un moyen pour spécifier les noeuds |
|    à partir desquels calculer la vraisemblance !!    |
\*----------------------------------------------------*/

void lik_dna_tree(int *n, int *s, int N, dna_matrix X, double *w,
		  matching match, double *edge_length, part_dna PART,
		  para_dna PARA, gamma_dna GAMMA, invar_dna INV,
		  double *BF, double *loglik)
{
    int i, j, k, nr, root;

    /* initialize before looping through the tree */
    /*-------------------------------------------------*\
    |     Cette boucle sera éventuellement ajustée      |
    | pour le calcul des vraisemblances conditionnelles |
    |      (sera même faite une fonction à part?)       |
    \*-------------------------------------------------*/
    nr = (*n*2 - 2);
    for(i = *n; i < *n*2 - 2; i++) { /* only for nodes */
        for(j = 0; j < *s; j++) {
	    k = i + j*nr;
	    X.A[k] = X.C[k] = X.G[k] = X.T[k] = 1.;
	}
    } /* end of initialization */

    /* loop through the matching */
    /*-------------------------------------------*\
    |  Ici aussi il faudra ajuster pour calculer  |
    | à partir des vraisemblances conditionnelles |
    \*-------------------------------------------*/
    /* We don't do the root node here,
       so k goes between 0 and n - 3. */
    for(k = 0; k < *n - 3; k++) {
        lik_dna_node(X, match.s1[k], match.s2[k],
		     match.a[k], edge_length, PART,
		     PARA, GAMMA, INV, BF, nr);
    }

    /* We now do the root */
    root = match.s2[*n - 2];
    lik_dna_root(X, match.s1[*n - 3], match.s2[*n - 3],
		 match.s1[*n - 2], root, edge_length,
                 PART, PARA, GAMMA, INV, BF, nr);
    *loglik = 0.;
    root--; /* to use root as index */
    for(j = 0; j < *s; j++) {
        i = root + j*nr;
        *loglik += w[j]*log(BF[0]*X.A[i] + BF[1]*X.C[i] + BF[2]*X.G[i] + BF[3]*X.T[i]);
    }
} /* lik_dna_tree */

double fcn_mlphylo_invar(double I, info *INFO)
{
    double loglik;
    int N;

    N = *(INFO->D->n)*2 - 3;
    INFO->D->INV.I[INFO->i] = I;
    lik_dna_tree(INFO->D->n, INFO->D->s, N, INFO->D->X, INFO->D->w,
		 INFO->D->match, INFO->D->edge_length,
		 INFO->D->partition, INFO->D->PAR, INFO->D->GAMMA,
		 INFO->D->INV, INFO->D->BF, &loglik);
    return(-loglik);
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
	D->INV.I[i] = I;
    }
}

double fcn_mlphylo_gamma(double a, info *INFO)
{
    double loglik;
    int N;

    N = *(INFO->D->n)*2 - 3;
    INFO->D->GAMMA.alpha[INFO->i] = a;

    lik_dna_tree(INFO->D->n, INFO->D->s, N, INFO->D->X, INFO->D->w,
		 INFO->D->match, INFO->D->edge_length,
		 INFO->D->partition, INFO->D->PAR, INFO->D->GAMMA,
		 INFO->D->INV, INFO->D->BF, &loglik);

    return(-loglik);
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
	D->GAMMA.alpha[i] = a;
    }
}

double fcn_mlphylo_para(double p, info *INFO)
{
    double loglik;
    int N;

    N = *(INFO->D->n)*2 - 3;
    INFO->D->PAR.para[INFO->i] = p;

    lik_dna_tree(INFO->D->n, INFO->D->s, N, INFO->D->X, INFO->D->w,
		 INFO->D->match, INFO->D->edge_length,
		 INFO->D->partition, INFO->D->PAR, INFO->D->GAMMA,
		 INFO->D->INV, INFO->D->BF, &loglik);

    return(-loglik);
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
	D->PAR.para[i] = p;
    }
}

double fcn_mlphylo_xi(double XI, info *INFO)
{
    double loglik;
    int N;

    N = *(INFO->D->n)*2 - 3;
    INFO->D->partition.xi[INFO->i] = XI;

    lik_dna_tree(INFO->D->n, INFO->D->s, N, INFO->D->X, INFO->D->w,
		 INFO->D->match, INFO->D->edge_length,
		 INFO->D->partition, INFO->D->PAR, INFO->D->GAMMA,
		 INFO->D->INV, INFO->D->BF, &loglik);

    return(-loglik);
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
	D->partition.xi[i] = XI;
    }
}

double fcn_mlphylo_edgelength(double l, info *INFO)
{
    double loglik;
    int N;

    N = *(INFO->D->n)*2 - 3;
    INFO->D->edge_length[INFO->i] = l;

    lik_dna_tree(INFO->D->n, INFO->D->s, N, INFO->D->X, INFO->D->w,
		 INFO->D->match, INFO->D->edge_length,
		 INFO->D->partition, INFO->D->PAR, INFO->D->GAMMA,
		 INFO->D->INV, INFO->D->BF, &loglik);
    return(-loglik);
}

void mlphylo_edgelength(int N, DNAdata *D, double *loglik)
/*
optimize branch lengths
*/
{
    int i;
    info INFO, *infptr;
    double l;

/*     int tmp; */

/*     tmp = D->GAMMA.ncat[0]; */
/*     D->GAMMA.ncat[0] = 1; */

    infptr = &INFO;
    INFO.D = D;

    for(i = 0; i < N; i++) {
        infptr->i = i;
	l = Brent_fmin(0.0, 0.1,
		       (double (*)(double, void*)) fcn_mlphylo_edgelength,
		       infptr, 1.e-6);
	D->edge_length[i] = l;
    }

/*     D->GAMMA.ncat[0] = tmp; */
}

/* void nni_matching(int *sib1, int *sib2, int *ancestor, */
/* 		  int node1, int node2, int crossed) */
/* { */
/*   int i, j, tmp; */

/*   /\* find the line with 'node1' as the ancestor *\/ */
/*   i = 0; */
/*   while (ancestor[i] != node1) i++; */
/*   /\* same thing for 'node2' *\/ */
/*   j = 0; */
/*   while (ancestor[j] != node1) j++; */

/*   /\* now do the inversion *\/ */
/*   if (crossed) { */
/*     tmp = sib2[i]; */
/*     sib2[i] = sib1[j]; */
/*     sib1[j] = tmp; */
/*   } else { */
/*     tmp = sib2[i]; */
/*     sib2[i] = sib2[j]; */
/*     sib2[j] = tmp; */
/*   } */

/*   /\* we check that both new pairs are correctly ordered *\/ */
/*   /\* (is this really needed?? or more checks??) *\/ */
/*   if (sib1[i] > sib2[i]) { */
/*     tmp = sib1[i]; */
/*     sib1[i] = sib2[i]; */
/*     sib2[i] = tmp; */
/*   } */
/*   if (sib1[j] > sib2[j]) { */
/*     tmp = sib1[j]; */
/*     sib1[j] = sib2[j]; */
/*     sib2[j] = tmp; */
/*   } */
/* } /\* nni_matching *\/ */

/* void mlphylo_topology(int n, DNAdata *D, double *loglik) */
/* { */
/*   int i, j, new_sib1, new_sib2, new_ancestor, *ns1, *ns2, *na; */
/*   double oldlik, newlik; */

/*   ns1 = &new_sib1; */
/*   ns2 = &new_sib2; */
/*   na = &new_ancestor; */

/*   ns1 = (int*)R_alloc(n - 1, sizeof(int)); */
/*   ns2 = (int*)R_alloc(n - 1, sizeof(int)); */
/*   na = (int*)R_alloc(n - 1, sizeof(int)); */

/*   oldlik = *loglik; */

/*   for (i = 0; i < n - 1; i++) { */
/*     ns1[i] = D->sib1[i]; */
/*     ns2[i] = D->sib2[i]; */
/*     na[i] = D->ancestor[i]; */
/*   } */

/*   for (i = 0; i < n - 3; i++) { */
/*     j = 0; */
/*     while (ns1[j] != na[i] || ns2[j] != na[i]) j++; */
/*     /\* peut-etre pas la pein de faire l'appel de fonction ou */
/*        alors, on peut le simplifier en passant juste le */
/*        numero de la ligne du matching *\/ */
/*     nni_matching(ns1, ns2, na, ns1[j], ns2[j], 0); */

/*     /\* In the followng, all the likelihood is re-computed: this *\/ */
/*     /\* may be speeded up by just re-computing the relevant part *\/ */
/*     /\* of the likelihood *\/ */
/*     lik_dna_tree(D->n, D->s, D->XA, D->XC, D->XG, D->XT, D->w, ns1, ns2, na, */
/* 		 D->edge_length, D->npart, D->partition, D->model, D->xi, D->para, */
/* 		 D->npara, D->pim_para, D->ncat, D->alpha, D->nalpha, D->pim_alpha, */
/* 		 D->invar, D->ninvar, D->pim_invar, D->BF, &newlik); */

/*     if (newlik < oldlik) { */
/*       for (i = 0; i < n - 1; i++) { */
/* 	ns1[i] = D->sib1[i]; */
/* 	ns2[i] = D->sib2[i]; */
/* 	na[i] = D->ancestor[i]; */
/*       } */
/*     } */
/*   } */
/* } /\* mlphylo_topology *\/ */

void mlphylo_DNAmodel(int *n, int *s, double *XA, double *XC, double *XG,
		     double *XT, double *w, int *sib1, int *sib2,
		     int *ancestor, double *edge_length, int *npart,
		     int *partition, int *model, double *xi, double *para,
		     int *npara, int *pim_para, double *alpha, int *nalpha,
		     int *pim_alpha, int *ncat, double *invar, int *ninvar,
		     int *pim_invar, double *BF, int *search_tree,
		     double *loglik)
/*
This function iterates to find the MLEs of the substitution
paramaters and of the branch lengths for a given tree.
*/
{
    int N, i, j, k;
    DNAdata *D, data;

    D = &data;
    N = *n*2 - 3;      /* the tree is unrooted */

    D->n = n;
    D->s = s;

    D->X.A = XA;
    D->X.C = XC;
    D->X.G = XG;
    D->X.T = XT;

    D->w = w;

    D->match.s1 = sib1;
    D->match.s2 = sib2;
    D->match.a = ancestor;

    D->edge_length = edge_length;

    D->partition.n = npart;
    D->partition.limit = partition;
    D->partition.model = model;
    D->partition.xi = xi;

    D->PAR.para = para;
    D->PAR.n = npara;
    D->PAR.pim = pim_para;

    D->GAMMA.ncat = ncat;
    D->GAMMA.alpha = alpha;
    D->GAMMA.n = nalpha;
    D->GAMMA.pim = pim_alpha;

    D->INV.I = invar;
    D->INV.n = ninvar;
    D->INV.pim = pim_invar;

    D->BF = BF;

    lik_dna_tree(D->n, D->s, N, D->X, D->w, D->match,
 		 D->edge_length, D->partition, D->PAR,
		 D->GAMMA, D->INV, D->BF, loglik);

/*     for (i = 0; i < 10; i++) { */
    if (*npart > 1) mlphylo_xi(*npart - 1, D, loglik);
    if (*npara) mlphylo_para(*npara, D, loglik);
    if (*nalpha) mlphylo_gamma(*nalpha, D, loglik);
    if (*ninvar) mlphylo_invar(*ninvar, D, loglik);
/*     mlphylo_edgelength(N, D, loglik); */

    lik_dna_tree(D->n, D->s, N, D->X, D->w, D->match,
 		 D->edge_length, D->partition, D->PAR,
		 D->GAMMA, D->INV, D->BF, loglik);

/*     if (*search_tree) mlphylo_topology(*n, D, loglik); */
/*     } */
} /* mlphylo_DNAmodel */

/*---------------------------------------------------------------*/

/* Exemples de "C-wrappers" pour obtenir la matrice */
/* de probabilités de transition. */

/* void titijc(double *t, double *u, double *P) */
/* { */
/*   PMAT_JC69(*t, *u, P); */
/* } */

/* void titihky(double *t, double *a, double *b, double *BF, double *P) */
/* { */
/*   PMAT_HKY85(*t, *a, *b, BF, P); */
/* } */
