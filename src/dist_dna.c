/* dist_dna.c       2008-12-22 */

/* Copyright 2005-2008 Emmanuel Paradis

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <R_ext/Lapack.h>

/* from R: print(log(4), d = 22) */
#define LN4 1.386294361119890572454

/* returns 8 if the base is known surely, 0 otherwise */
#define KnownBase(a) a & 8

/* returns 1 if the base is adenine surely, 0 otherwise */
#define IsAdenine(a) a == 136

/* returns 1 if the base is guanine surely, 0 otherwise */
#define IsGuanine(a) a == 72

/* returns 1 if the base is cytosine surely, 0 otherwise */
#define IsCytosine(a) a == 40

/* returns 1 if the base is thymine surely, 0 otherwise */
#define IsThymine(a) a == 24

/* returns 1 if the base is a purine surely, 0 otherwise */
#define IsPurine(a) a > 63

/* returns 1 if the base is a pyrimidine surely, 0 otherwise */
#define IsPyrimidine(a) a < 64

/* returns 1 if both bases are different surely, 0 otherwise */
#define DifferentBase(a, b) (a & b) < 16

/* returns 1 if both bases are the same surely, 0 otherwise */
#define SameBase(a, b) KnownBase(a) && a == b

/* computes directly the determinant of a 4x4 matrix */
double detFourByFour(double *x)
{
    double det, a33a44, a34a43, a34a42, a32a44, a32a43, a33a42, a34a41, a31a44, a31a43, a33a41, a31a42, a32a41;

    a33a44 = x[10]*x[15]; a34a43 = x[14]*x[11];
    a34a42 = x[14]*x[7];  a32a44 = x[6]*x[15];
    a32a43 = x[6]*x[11];  a33a42 = x[10]*x[7];
    a34a41 = x[14]*x[3];  a31a44 = x[2]*x[15];
    a31a43 = x[2]*x[11];  a33a41 = x[10]*x[3];
    a31a42 = x[2]*x[7];   a32a41 = x[6]*x[3];

    det = x[0]*x[5]*(a33a44 - a34a43) + x[0]*x[9]*(a34a42 - a32a44) +
      x[0]*x[13]*(a32a43 - a33a42) + x[4]*x[9]*(a31a44 - a34a41) +
      x[4]*x[13]*(a33a41 - a31a43) + x[4]*x[1]*(a34a43 - a33a44) +
      x[8]*x[13]*(a31a42 - a32a41) + x[8]*x[1]*(a32a44 - a34a42) +
      x[8]*x[5]*(a34a41 - a31a44) + x[12]*x[1]*(a33a42 - a32a43) +
      x[12]*x[5]*(a31a43 - a33a41) + x[12]*x[9]*(a32a41 - a31a42);

    return det;
}

#define CHECK_PAIRWISE_DELETION\
    if (KnownBase(x[s1]) && KnownBase(x[s2])) L++;\
    else continue;

void distDNA_raw(unsigned char *x, int *n, int *s, double *d, int scaled)
{
    int i1, i2, s1, s2, target, Nd;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n)
	      if (DifferentBase(x[s1], x[s2])) Nd++;
	    if (scaled) d[target] = ((double) Nd / *s);
	    else d[target] = ((double) Nd);
	    target++;
	}
    }
}

void distDNA_raw_pairdel(unsigned char *x, int *n, int *s, double *d, int scaled)
{
    int i1, i2, s1, s2, target, Nd, L;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
                CHECK_PAIRWISE_DELETION
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    if (scaled) d[target] = ((double) Nd/L);
	    else d[target] = ((double) Nd);
	    target++;
	}
    }
}

#define COMPUTE_DIST_JC69\
    p = ((double) Nd/L);\
    if (*gamma)\
      d[target] = 0.75 * *alpha*(pow(1 - 4*p/3, -1/ *alpha) - 1);\
    else d[target] = -0.75*log(1 - 4*p/3);\
    if (*variance) {\
        if (*gamma) var[target] = p*(1 - p)/(pow(1 - 4*p/3, -2/(*alpha + 1)) * L);\
	else var[target] = p*(1 - p)/(pow(1 - 4*p/3, 2)*L);\
    }

void distDNA_JC69(unsigned char *x, int *n, int *s, double *d,
		  int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p;

    L = *s;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n)
	      if (DifferentBase(x[s1], x[s2])) Nd++;
	    COMPUTE_DIST_JC69
	    target++;
	}
    }
}

void distDNA_JC69_pairdel(unsigned char *x, int *n, int *s, double *d,
			  int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    COMPUTE_DIST_JC69
	    target++;
	}
    }
}

#define COUNT_TS_TV\
    if (SameBase(x[s1], x[s2])) continue;\
    Nd++;\
    if (IsPurine(x[s1]) && IsPurine(x[s2])) {\
        Ns++;\
        continue;\
    }\
    if (IsPyrimidine(x[s1]) && IsPyrimidine(x[s2])) Ns++;

#define COMPUTE_DIST_K80\
    P = ((double) Ns/L);\
    Q = ((double) (Nd - Ns)/L);\
    a1 = 1 - 2*P - Q;\
    a2 = 1 - 2*Q;\
    if (*gamma) {\
        b = -1 / *alpha;\
    	d[target] = *alpha * (pow(a1, b) + 0.5*pow(a2, b) - 1.5)/2;\
    }\
    else d[target] = -0.5 * log(a1 * sqrt(a2));\
    if (*variance) {\
        if (*gamma) {\
    	    b = -(1 / *alpha + 1);\
    	    c1 = pow(a1, b);\
    	    c2 = pow(a2, b);\
    	    c3 = (c1 + c2)/2;\
    	} else {\
    	  c1 = 1/a1;\
    	  c2 = 1/a2;\
    	  c3 = (c1 + c2)/2;\
    	}\
    	var[target] = (c1*c1*P + c3*c3*Q - pow(c1*P + c3*Q, 2))/L;\
    }

void distDNA_K80(unsigned char *x, int *n, int *s, double *d,
		 int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, Ns, L;
    double P, Q, a1, a2, b, c1, c2, c3;

    L = *s;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_K80
	    target++;
	}
    }
}

void distDNA_K80_pairdel(unsigned char *x, int *n, int *s, double *d,
			 int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, Ns, L;
    double P, Q, a1, a2, b, c1, c2, c3;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
		COUNT_TS_TV
	    }
	    COMPUTE_DIST_K80
	    target++;
	}
    }
}

#define COMPUTE_DIST_F81\
    p = ((double) Nd/L);\
    if (*gamma) d[target] = E * *alpha * (pow(1 - p/E, -1/ *alpha) - 1);\
    else d[target] = -E*log(1 - p/E);\
    if (*variance) {\
	if (*gamma) var[target] = p*(1 - p)/(pow(1 - p/E, -2/(*alpha + 1)) * L);\
	else var[target] = p*(1 - p)/(pow(1 - p/E, 2)*L);\
    }

void distDNA_F81(unsigned char *x, int *n, int *s, double *d, double *BF,
		 int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p, E;

    L = *s;
    E = 1 - BF[0]*BF[0] - BF[1]*BF[1] - BF[2]*BF[2] - BF[3]*BF[3];

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n)
	      if (DifferentBase(x[s1], x[s2])) Nd++;
	    COMPUTE_DIST_F81
	    target++;
	}
    }
}

void distDNA_F81_pairdel(unsigned char *x, int *n, int *s, double *d, double *BF,
			 int *variance, double *var, int *gamma, double *alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p, E;

    E = 1 - BF[0]*BF[0] - BF[1]*BF[1] - BF[2]*BF[2] - BF[3]*BF[3];

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    COMPUTE_DIST_F81
	    target++;
	}
    }
}

#define COUNT_TS_TV1_TV2\
    if (SameBase(x[s1], x[s2])) continue;\
    Nd++;\
    if ((x[s1] | x[s2]) == 152 || (x[s1] | x[s2]) == 104) {\
        Nv1++;\
        continue;\
    }\
    if ((x[s1] | x[s2]) == 168 || (x[s1] | x[s2]) == 88) Nv2++;


#define COMPUTE_DIST_K81\
    P = ((double) (Nd - Nv1 - Nv2)/L);\
    Q = ((double) Nv1/L);\
    R = ((double) Nv2/L);\
    a1 = 1 - 2*P - 2*Q;\
    a2 = 1 - 2*P - 2*R;\
    a3 = 1 - 2*Q - 2*R;\
    d[target] = -0.25*log(a1*a2*a3);\
    if (*variance) {\
        a = (1/a1 + 1/a2)/2;\
    	b = (1/a1 + 1/a3)/2;\
    	c = (1/a2 + 1/a3)/2;\
      var[target] = (a*a*P + b*b*Q + c*c*R - pow(a*P + b*Q + c*R, 2))/2;\
    }

void distDNA_K81(unsigned char *x, int *n, int *s, double *d,
		 int *variance, double *var)
{
    int i1, i2, Nd, Nv1, Nv2, L, s1, s2, target;
    double P, Q, R, a1, a2, a3, a, b, c;

    L = *s;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = Nv1 = Nv2 = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        COUNT_TS_TV1_TV2
	    }
	    COMPUTE_DIST_K81
	    target++;
	}
    }
}

void distDNA_K81_pairdel(unsigned char *x, int *n, int *s, double *d,
			 int *variance, double *var)
{
    int i1, i2, Nd, Nv1, Nv2, L, s1, s2, target;
    double P, Q, R, a1, a2, a3, a, b, c;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
  	    Nd = Nv1 = Nv2 = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
	        COUNT_TS_TV1_TV2
	    }
	    COMPUTE_DIST_K81
	    target++;
	}
    }
}

#define PREPARE_BF_F84\
    A = (BF[0]*BF[2])/(BF[0] + BF[2]) + (BF[1]*BF[3])/(BF[1] + BF[3]);\
    B = BF[0]*BF[2] + BF[1]*BF[3];\
    C = (BF[0] + BF[2])*(BF[1] + BF[3]);

#define COMPUTE_DIST_F84\
   P = ((double) Ns/L);\
   Q = ((double) (Nd - Ns)/L);\
   d[target] = -2*A*log(1 - (P/(2*A) - (A - B)*Q/(2*A*C))) + 2*(A - B - C)*log(1 - Q/(2*C));\
   if (*variance) {\
       t1 = A*C;\
       t2 = C*P/2;\
       t3 = (A - B)*Q/2;\
       a = t1/(t1 - t2 - t3);\
       b = A*(A - B)/(t1 - t2 - t3) - (A - B - C)/(C - Q/2);\
       var[target] = (a*a*P + b*b*Q - pow(a*P + b*Q, 2))/2;\
   }

void distDNA_F84(unsigned char *x, int *n, int *s, double *d,
		 double *BF, int *variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, A, B, C, a, b, t1, t2, t3;

    PREPARE_BF_F84
    L = *s;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_F84
	    target++;
	}
    }
}

void distDNA_F84_pairdel(unsigned char *x, int *n, int *s, double *d,
			 double *BF, int *variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, A, B, C, a, b, t1, t2, t3;

    PREPARE_BF_F84

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		CHECK_PAIRWISE_DELETION
		COUNT_TS_TV
	    }
	    COMPUTE_DIST_F84
	    target++;
	}
    }
}

#define COMPUTE_DIST_T92\
    P = ((double) Ns/L);\
    Q = ((double) (Nd - Ns)/L);\
    a1 = 1 - P/wg - Q;\
    a2 = 1 - 2*Q;\
    d[target] = -wg*log(a1) - 0.5*(1 - wg)*log(a2);\
    if (*variance) {\
        c1 = 1/a1;\
        c2 = 1/a2;\
        c3 = wg*(c1 - c2) + c2;\
        var[target] = (c1*c1*P + c3*c3*Q - pow(c1*P + c3*Q, 2))/L;\
    }

void distDNA_T92(unsigned char *x, int *n, int *s, double *d,
		 double *BF, int *variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, wg, a1, a2, c1, c2, c3;

    L = *s;
    wg = 2 * (BF[1] + BF[2]) * (1 - (BF[1] + BF[2]));

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_T92
	    target++;
	}
    }
}

void distDNA_T92_pairdel(unsigned char *x, int *n, int *s, double *d,
			 double *BF, int *variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, wg, a1, a2, c1, c2, c3;

    wg = 2 * (BF[1] + BF[2]) * (1 - (BF[1] + BF[2]));

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_T92
	    target++;
	}
    }
}

/* returns 1 if one of the base is adenine and
   the other one is guanine surely, 0 otherwise */
#define AdenineAndGuanine(a, b) (a | b) == 200

/* returns 1 if one of the base is cytosine and
   the other one is thymine surely, 0 otherwise */
#define CytosineAndThymine(a, b) (a | b) == 56

#define PREPARE_BF_TN93\
    gR = BF[0] + BF[2];\
    gY = BF[1] + BF[3];\
    k1 = 2 * BF[0] * BF[2] / gR;\
    k2 = 2 * BF[1] * BF[3] / gY;\
    k3 = 2 * (gR * gY - BF[0]*BF[2]*gY/gR - BF[1]*BF[3]*gR/gY);

#define COUNT_TS1_TS2_TV\
    if (DifferentBase(x[s1], x[s2])) {\
        Nd++;\
        if (AdenineAndGuanine(x[s1], x[s2])) {\
            Ns1++;\
    	continue;\
        }\
        if (CytosineAndThymine(x[s1], x[s2])) Ns2++;\
    }

#define COMPUTE_DIST_TN93\
    P1 = ((double) Ns1/L);\
    P2 = ((double) Ns2/L);\
    Q = ((double) (Nd - Ns1 - Ns2)/L);\
    w1 = 1 - P1/k1 - Q/(2*gR);\
    w2 = 1 - P2/k2 - Q/(2*gY);\
    w3 = 1 - Q/(2*gR*gY);\
    if (*gamma) {\
        k4 = 2*(BF[0]*BF[2] + BF[1]*BF[3] + gR*gY);\
    	b = -1 / *alpha;\
    	c1 = pow(w1, b);\
    	c2 = pow(w2, b);\
    	c3 = pow(w3, b);\
    	c4 = k1*c1/(2*gR) + k2*c2/(2*gY) + k3*c3/(2*gR*gY);\
    	d[target] = *alpha * (k1*pow(w1, b) + k2*pow(w2, b) + k3*pow(w3, b) - k4);\
    } else {\
        k4 = 2*((BF[0]*BF[0] + BF[2]*BF[2])/(2*gR*gR) + (BF[2]*BF[2] + BF[3]*BF[3])/(2*gY*gY));\
    	c1 = 1/w1;\
    	c2 = 1/w2;\
    	c3 = 1/w3;\
    	c4 = k1 * c1/(2 * gR) + k2 * c2/(2 * gY) + k4 * c3;\
    	d[target] = -k1*log(w1) - k2*log(w2) - k3*log(w3);\
    }\
    if (*variance)\
      var[target] = (c1*c1*P1 + c2*c2*P2 + c4*c4*Q - pow(c1*P1 + c2*P2 + c4*Q, 2))/L;

void distDNA_TN93(unsigned char *x, int *n, int *s, double *d,
		  double *BF, int *variance, double *var,
		  int *gamma, double *alpha)
{
    int i1, i2, k, Nd, Ns1, Ns2, L, target, s1, s2;
    double P1, P2, Q, A, B, C, gR, gY, k1, k2, k3, k4, w1, w2, w3, c1, c2, c3, c4, b;

    L = *s;

    PREPARE_BF_TN93

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns1 = Ns2 = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		COUNT_TS1_TS2_TV
	    }
	    COMPUTE_DIST_TN93
	    target++;
	}
    }
}

void distDNA_TN93_pairdel(unsigned char *x, int *n, int *s, double *d,
			  double *BF, int *variance, double *var,
			  int *gamma, double *alpha)
{
    int i1, i2, k, Nd, Ns1, Ns2, L, target, s1, s2;
    double P1, P2, Q, A, B, C, gR, gY, k1, k2, k3, k4, w1, w2, w3, c1, c2, c3, c4, b;

    PREPARE_BF_TN93

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns1 = Ns2 = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		CHECK_PAIRWISE_DELETION
		COUNT_TS1_TS2_TV
	    }
	    COMPUTE_DIST_TN93
	    target++;
	}
    }
}

void distDNA_GG95(unsigned char *x, int *n, int *s, double *d,
		  int *variance, double *var)
{
    int i1, i2, s1, s2, target, GC, Nd, Ns, tl, npair;
    double *theta, gcprop, *P, pp, *Q, qq, *tstvr, svr, A, sum, ma /* mean alpha */, K1, K2;

    theta = &gcprop;
    P = &pp;
    Q = &qq;
    tstvr = &svr;

    npair = *n * (*n - 1) / 2;

    theta = (double*)R_alloc(*n, sizeof(double));
    P = (double*)R_alloc(npair, sizeof(double));
    Q = (double*)R_alloc(npair, sizeof(double));
    tstvr = (double*)R_alloc(npair, sizeof(double));

    /* get the proportion of GC (= theta) in each sequence */
    for (i1 = 1; i1 <= *n; i1++) {
        GC = 0;
	for (s1 = i1 - 1; s1 < i1 + *n*(*s - 1); s1 += *n)
	  if (IsCytosine(x[s1]) || IsGuanine(x[s1])) GC += 1;
	theta[i1 - 1] = ((double) GC / *s);
    }

    /* get the proportions of transitions and transversions,
       and the estimates of their ratio for each pair */
    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        COUNT_TS_TV
	    }
	    P[target] = ((double) Ns / *s);
	    Q[target] = ((double) (Nd - Ns) / *s);
	    A = log(1 - 2*Q[target]);
	    tstvr[target] = 2*(log(1 - 2*P[target] - Q[target]) - 0.5*A)/A;
	    target++;
	}
    }

    /* compute the mean alpha (ma) = mean Ts/Tv */
    sum = 0;
    tl = 0;
    for (i1 = 0; i1 < npair; i1++)
    /* some values of tstvr are -Inf if there is no
       transversions observed: we exclude them */
      if (R_FINITE(tstvr[i1])) {
	  sum += tstvr[i1];
	  tl += 1;
      }
    ma = sum/tl;

    /* compute the distance for each pair */
    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    A = 1 - 2*Q[target];
	    K1 = 1 + ma*(theta[i1 - 1]*(1 - theta[i1 - 1]) + theta[i2 - 1]*(1 - theta[i2 - 1]));
	    K2 = ma*pow(theta[i1 - 1] - theta[i2 - 1], 2)/(ma + 1);
	    d[target] = -0.5*K1*log(A) + K2*(1 - pow(A, 0.25*(ma + 1)));
	    if (*variance)
	      var[target] = pow(K1 + K2*0.5*(ma + 1)*pow(A, 0.25*(ma + 1)), 2)*Q[target]*(1 - Q[target])/(A*A * *s);
	    target++;
	}
    }
}

void distDNA_GG95_pairdel(unsigned char *x, int *n, int *s, double *d,
			  int *variance, double *var)
{
    int i1, i2, s1, s2, target, *L, length, GC, Nd, Ns, tl, npair;
    double *theta, gcprop, *P, pp, *Q, qq, *tstvr, svr, A, sum, ma /* mean alpha */, K1, K2;

    theta = &gcprop;
    L = &length;
    P = &pp;
    Q = &qq;
    tstvr = &svr;

    npair = *n * (*n - 1) / 2;

    theta = (double*)R_alloc(*n, sizeof(double));
    L = (int*)R_alloc(npair, sizeof(int));
    P = (double*)R_alloc(npair, sizeof(double));
    Q = (double*)R_alloc(npair, sizeof(double));
    tstvr = (double*)R_alloc(npair, sizeof(double));

    /* get the proportion of GC (= theta) in each sequence */
    for (i1 = 1; i1 <= *n; i1++) {
        tl = GC = 0;
	for (s1 = i1 - 1; s1 < i1 + *n*(*s - 1); s1 += *n) {
	    if (KnownBase(x[s1])) tl++;
	    else continue;
	    if (IsCytosine(x[s1]) || IsGuanine(x[s1])) GC += 1;
	}
	theta[i1 - 1] = ((double) GC / tl);
    }

    /* get the proportions of transitions and transversions,
       and the estimates of their ratio for each pair; we
       also get the sample size for each pair in L */
    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    Nd = Ns = L[target] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        if (KnownBase(x[s1]) && KnownBase(x[s2])) L[target]++;
		else continue;
	        COUNT_TS_TV
	    }
	    P[target] = ((double) Ns/L[target]);
	    Q[target] = ((double) (Nd - Ns)/L[target]);
	    A = log(1 - 2*Q[target]);
	    tstvr[target] = 2*(log(1 - 2*P[target] - Q[target]) - 0.5*A)/A;
	    target++;
	}
    }

    /* compute the mean alpha (ma) = mean Ts/Tv */
    sum = 0;
    tl = 0;
    for (i1 = 0; i1 < npair; i1++)
    /* some values of tstvr are -Inf if there is no
       transversions observed: we exclude them */
      if (R_FINITE(tstvr[i1])) {
	  sum += tstvr[i1];
	  tl += 1;
      }
    ma = sum/tl;

    /* compute the distance for each pair */
    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    A = 1 - 2*Q[target];
	    K1 = 1 + ma*(theta[i1 - 1]*(1 - theta[i1 - 1]) + theta[i2 - 1]*(1 - theta[i2 - 1]));
	    K2 = ma*pow(theta[i1 - 1] - theta[i2 - 1], 2)/(ma + 1);
	    d[target] = -0.5*K1*log(A) + K2*(1 - pow(A, 0.25*(ma + 1)));
	    if (*variance)
	      var[target] = pow(K1 + K2*0.5*(ma + 1)*pow(A, 0.25*(ma + 1)), 2)*Q[target]*(1 - Q[target])/(A*A*L[target]);
	    target++;
	}
    }
}

#define DO_CONTINGENCY_NUCLEOTIDES\
    switch (x[s1]) {\
    case 136 : m = 0; break;\
    case 72 : m = 1; break;\
    case 40 : m = 2; break;\
    case 24 : m = 3; break;\
    }\
    switch (x[s2]) {\
    case 72 : m += 4; break;\
    case 40 : m += 8; break;\
    case 24 : m += 12; break;\
    }\
    Ntab[m]++;

#define COMPUTE_DIST_LogDet\
    for (k = 0; k < 16; k++) Ftab[k] = ((double) Ntab[k]/L);\
    d[target] = (-log(detFourByFour(Ftab))/4 - LN4)/4;\
    if (*variance) {\
        /* For the inversion, we first make U an identity matrix */\
        for (k = 1; k < 15; k++) U[k] = 0;\
    	U[0] = U[5] = U[10] = U[15] = 1;\
    	/* The matrix is not symmetric, so we use 'dgesv'. */\
    	/* This subroutine puts the result in U. */\
    	F77_CALL(dgesv)(&ndim, &ndim, Ftab, &ndim, ipiv, U, &ndim, &info);\
    	var[target] = (U[0]*U[0]*Ftab[0] + U[1]*U[1]*Ftab[4] +\
    		       U[2]*U[2]*Ftab[8] + U[3]*U[3]*Ftab[12] +\
    		       U[4]*U[4]*Ftab[1] + U[5]*U[5]*Ftab[5] +\
    		       U[6]*U[6]*Ftab[9] + U[7]*U[7]*Ftab[13] +\
    		       U[8]*U[8]*Ftab[2] + U[9]*U[9]*Ftab[6] +\
    		       U[10]*U[10]*Ftab[10] + U[11]*U[11]*Ftab[14] +\
    		       U[12]*U[12]*Ftab[3] + U[13]*U[13]*Ftab[7] +\
    		       U[14]*U[14]*Ftab[11] + U[15]*U[15]*Ftab[15] - 16)/(L*16);\
    }

void distDNA_LogDet(unsigned char *x, int *n, int *s, double *d,
		    int *variance, double *var)
{
    int i1, i2, k, m, s1, s2, target, L, Ntab[16], ndim = 4, info, ipiv[16];
    double Ftab[16], U[16];

    L = *s;

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		DO_CONTINGENCY_NUCLEOTIDES
	    }
	    COMPUTE_DIST_LogDet
	    target++;
	}
    }
}

void distDNA_LogDet_pairdel(unsigned char *x, int *n, int *s, double *d,
			    int *variance, double *var)
{
    int i1, i2, k, m, s1, s2, target, L, Ntab[16], ndim = 4, info, ipiv[16];
    double Ftab[16], U[16];

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		CHECK_PAIRWISE_DELETION
		DO_CONTINGENCY_NUCLEOTIDES
	    }
	    COMPUTE_DIST_LogDet
	    target++;
	}
    }
}

void distDNA_BH87(unsigned char *x, int *n, int *s, double *d,
		  int *variance, double *var)
/* <FIXME>
   For the moment there is no need to check for pairwise deletions
   since DO_CONTINGENCY_NUCLEOTIDES considers only the known nucleotides.
   In effect the pairwise deletion has possibly been done before.
   The sequence length(s) are used only to compute the variances, which is
   currently not available.
   </FIXME> */
{
    int i1, i2, k, kb, s1, s2, m, Ntab[16], ROWsums[4], ndim = 4, info, ipiv[16];
    double P12[16], P21[16], U[16];

    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		DO_CONTINGENCY_NUCLEOTIDES
	    }

            /* get the rowwise sums of Ntab */
            ROWsums[0] = Ntab[0] + Ntab[4] + Ntab[8] + Ntab[12];
            ROWsums[1] = Ntab[1] + Ntab[5] + Ntab[9] + Ntab[13];
            ROWsums[2] = Ntab[2] + Ntab[6] + Ntab[10] + Ntab[14];
            ROWsums[3] = Ntab[3] + Ntab[7] + Ntab[11] + Ntab[15];

            for (k = 0; k < 16; k++)
              P12[k] = ((double) Ntab[k]);

            /* scale each element of P12 by its rowwise sum */
            for (k = 0; k < 4; k++)
              for (kb = 0; kb < 16; kb += 4)
            	P12[k + kb] = P12[k + kb]/ROWsums[k];

            d[*n*(i2 - 1) + i1 - 1] = -log(detFourByFour(P12))/4;

            /* compute the columnwise sums of Ntab: these
               are the rowwise sums of its transpose */
            ROWsums[0] = Ntab[0] + Ntab[1] + Ntab[2] + Ntab[3];
            ROWsums[1] = Ntab[4] + Ntab[5] + Ntab[6] + Ntab[7];
            ROWsums[2] = Ntab[8] + Ntab[9] + Ntab[10] + Ntab[11];
            ROWsums[3] = Ntab[12] + Ntab[13] + Ntab[14] + Ntab[15];

            /* transpose Ntab and store the result in P21 */
            for (k = 0; k < 4; k++)
               for (kb = 0; kb < 4; kb++)
            	 P21[kb + 4*k] = Ntab[k + 4*kb];

            /* scale as above */
            for (k = 0; k < 4; k++)
              for (kb = 0; kb < 16; kb += 4)
            	P21[k + kb] = P21[k + kb]/ROWsums[k];

            d[*n*(i1 - 1) + i2 - 1] = -log(detFourByFour(P21))/4;
	}
    }
}

#define COMPUTE_DIST_ParaLin\
    for (k = 0; k < 16; k++) Ftab[k] = ((double) Ntab[k]/L);\
    d[target] = -log(detFourByFour(Ftab)/\
		     sqrt(find[0][i1 - 1]*find[1][i1 - 1]*find[2][i1 - 1]*find[3][i1 - 1]*\
			  find[0][i2 - 1]*find[1][i2 - 1]*find[2][i2 - 1]*find[3][i2 - 1]))/4;\
    if (*variance) {\
        /* For the inversion, we first make U an identity matrix */\
        for (k = 1; k < 15; k++) U[k] = 0;\
    	U[0] = U[5] = U[10] = U[15] = 1;\
    	/* The matrix is not symmetric, so we use 'dgesv'. */\
    	/* This subroutine puts the result in U. */\
    	F77_CALL(dgesv)(&ndim, &ndim, Ftab, &ndim, ipiv, U, &ndim, &info);\
    	var[target] = (U[0]*U[0]*Ftab[0] + U[1]*U[1]*Ftab[4] +\
    		       U[2]*U[2]*Ftab[8] + U[3]*U[3]*Ftab[12] +\
    		       U[4]*U[4]*Ftab[1] + U[5]*U[5]*Ftab[5] +\
    		       U[6]*U[6]*Ftab[9] + U[7]*U[7]*Ftab[13] +\
    		       U[8]*U[8]*Ftab[2] + U[9]*U[9]*Ftab[6] +\
    		       U[10]*U[10]*Ftab[10] + U[11]*U[11]*Ftab[14] +\
    		       U[12]*U[12]*Ftab[3] + U[13]*U[13]*Ftab[7] +\
    		       U[14]*U[14]*Ftab[11] + U[15]*U[15]*Ftab[15] -\
		       4*(1/sqrt(find[0][i1 - 1]*find[0][i2 - 1]) +\
                       1/sqrt(find[1][i1 - 1]*find[1][i2 - 1]) +\
		       1/sqrt(find[2][i1 - 1]*find[2][i2 - 1]) +\
                       1/sqrt(find[3][i1 - 1]*find[3][i2 - 1])))/(L*16);\
    }

void distDNA_ParaLin(unsigned char *x, int *n, int *s, double *d,
		     int *variance, double *var)
{
    int i1, i2, k, s1, s2, m, target, L, Ntab[16], ndim = 4, info, ipiv[16];
    double Ftab[16], U[16], *find[4];

    L = *s;

    for (k = 0; k < 4; k++)
      find[k] = (double*)R_alloc(*n, sizeof(double));

    for (i1 = 0; i1 < *n; i1++)
      for (k = 0; k < 4; k++) find[k][i1] = 0.0;

    for (i1 = 0; i1 < *n; i1++) {
        for (s1 = i1; s1 < i1 + *n*(*s - 1) + 1; s1+= *n) {
            switch (x[s1]) {
	    case 136 : find[0][i1]++; break;
	    case 40 : find[1][i1]++; break;
	    case 72 : find[2][i1]++; break;
	    case 24 : find[3][i1]++; break;
	    }
        }
        for (k = 0; k < 4; k++) find[k][i1] /= L;
    }

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
		DO_CONTINGENCY_NUCLEOTIDES
	    }
	    COMPUTE_DIST_ParaLin
	    target++;
	}
    }
}

void distDNA_ParaLin_pairdel(unsigned char *x, int *n, int *s, double *d,
			     int *variance, double *var)
{
    int i1, i2, k, s1, s2, m, target, L, Ntab[16], ndim = 4, info, ipiv[16];
    double Ftab[16], U[16], *find[4];

    L = 0;

    for (k = 0; k < 4; k++)
      find[k] = (double*)R_alloc(*n, sizeof(double));

    for (i1 = 0; i1 < *n; i1++)
      for (k = 0; k < 4; k++) find[k][i1] = 0.0;

    for (i1 = 0; i1 < *n; i1++) {
        L = 0;
        for (s1 = i1; s1 < i1 + *n*(*s - 1) + 1; s1+= *n) {
	    if (KnownBase(x[s1])) {
	        L++;
                switch (x[s1]) {
	        case 136 : find[0][i1]++; break;
	        case 40 : find[1][i1]++; break;
	        case 72 : find[2][i1]++; break;
	        case 24 : find[3][i1]++; break;
	        }
	    }
        }
        for (k = 0; k < 4; k++) find[k][i1] /= L;
    }

    target = 0;
    for (i1 = 1; i1 < *n; i1++) {
        for (i2 = i1 + 1; i2 <= *n; i2++) {
	    L = 0;
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + *n*(*s - 1); s1+= *n, s2 += *n) {
	        CHECK_PAIRWISE_DELETION
		DO_CONTINGENCY_NUCLEOTIDES
	    }
	    COMPUTE_DIST_ParaLin
	    target++;
	}
    }
}

void BaseProportion(unsigned char *x, int *n, double *BF)
{
    int i, m;

    m = 0;
    for (i = 0; i < *n; i++) {
        if (KnownBase(x[i])) {
	    m++;
	    switch (x[i]) {
	    case 136 : BF[0]++; break;
	    case 40 : BF[1]++; break;
	    case 72 : BF[2]++; break;
	    case 24 : BF[3]++; break;
	    }
	}
    }
    for (i = 0; i < 4; i++) BF[i] /= m;
}

void SegSites(unsigned char *x, int *n, int *s, int *seg)
{
    int i, j;
    unsigned char basis;

    for (j = 0; j < *s; j++) {
        i = *n * j;
	while (!KnownBase(x[i])) i++;
	basis = x[i];
	i++;
	while (i < *n * (j + 1)) {
	    if (x[i] == basis) i++;
	    else {
	        seg[j] = 1;
		break;
	    }
	}
    }
}

void GlobalDeletionDNA(unsigned char *x, int *n, int *s, int *keep)
{
    int i, j;

    for (j = 0; j < *s; j++) {
        i = *n * j;
	while (i < *n * (j + 1)) {
	    if (KnownBase(x[i])) i++;
	    else {
	        keep[j] = 0;
		break;
	    }
	}
    }
}

void dist_dna(unsigned char *x, int *n, int *s, int *model, double *d,
	      double *BF, int *pairdel, int *variance, double *var,
	      int *gamma, double *alpha)
{
    switch (*model) {
    case 1 : if (pairdel) distDNA_raw_pairdel(x, n, s, d, 1);
             else distDNA_raw(x, n, s, d, 1); break;

    case 2 : if (pairdel) distDNA_JC69_pairdel(x, n, s, d, variance, var, gamma, alpha);
             else distDNA_JC69(x, n, s, d, variance, var, gamma, alpha); break;

    case 3 : if (pairdel) distDNA_K80_pairdel(x, n, s, d, variance, var, gamma, alpha);
             else distDNA_K80(x, n, s, d, variance, var, gamma, alpha); break;

    case 4 : if (pairdel) distDNA_F81_pairdel(x, n, s, d, BF, variance, var, gamma, alpha);
             else distDNA_F81(x, n, s, d, BF, variance, var, gamma, alpha); break;

    case 5 : if (pairdel) distDNA_K81_pairdel(x, n, s, d, variance, var);
             else distDNA_K81(x, n, s, d, variance, var); break;

    case 6 : if (pairdel) distDNA_F84_pairdel(x, n, s, d, BF, variance, var);
             else distDNA_F84(x, n, s, d, BF, variance, var); break;

    case 7 : if (pairdel) distDNA_T92_pairdel(x, n, s, d, BF, variance, var);
             else distDNA_T92(x, n, s, d, BF, variance, var); break;

    case 8 : if (pairdel) distDNA_TN93_pairdel(x, n, s, d, BF, variance, var, gamma, alpha);
             else distDNA_TN93(x, n, s, d, BF, variance, var, gamma, alpha); break;

    case 9 : if (pairdel) distDNA_GG95_pairdel(x, n, s, d, variance, var);
             else distDNA_GG95(x, n, s, d, variance, var); break;

    case 10 : if (pairdel) distDNA_LogDet_pairdel(x, n, s, d, variance, var);
              else distDNA_LogDet(x, n, s, d, variance, var); break;

    case 11 : distDNA_BH87(x, n, s, d, variance, var); break;

    case 12 : if (pairdel) distDNA_ParaLin_pairdel(x, n, s, d, variance, var);
              else distDNA_ParaLin(x, n, s, d, variance, var); break;
    case 13 : if (pairdel) distDNA_raw_pairdel(x, n, s, d, 0);
             else distDNA_raw(x, n, s, d, 0); break;
    }
}
