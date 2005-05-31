/* Copyright 2005 Emmanuel Paradis <paradis@isem.univ-montp2.fr> */

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

#define both_non_n(a, b) (strcmp(a, "n") && strcmp(b, "n"))

void dist_dna_JC69(char **x, int *n, int *s, double *d, int *pairdel,
		   int *variance, double *var, int *gamma, double *alpha)
{
  int i, j, k, s1, s2, target, Nd, L;
  double p;

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) Nd += 1;
      }
      if (!*pairdel) L = *s;
      p = ((double) Nd/L);
      if (*gamma) d[target] = 0.75 * *alpha * (pow(1 - 4*p/3, -1/ *alpha) - 1);
      else d[target] = -0.75 * log(1 - 4 * p/3);
      if (*variance) {
	if (*gamma) var[target] = p*(1 - p)/(pow(1 - 4*p/3, -2/(*alpha + 1)) * L);
	else var[target] = p*(1 - p)/(pow(1 - 4*p/3, 2)*L);
      }
    }
  }
}

void dist_dna_K80(char **x, int *n, int *s, double *d, int *pairdel,
		   int *variance, double *var, int *gamma, double *alpha)
{
  int i, j, k, s1, s2, target, Nd, Ns, L;
  double P, Q, a1, a2, b, c1, c2, c3;

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Ns = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a") && !strcmp(x[s2], "g"))
	    Ns += 1;
	  else {
	    if (!strcmp(x[s1], "g") && !strcmp(x[s2], "a"))
	      Ns += 1;
	    else {
	      if (!strcmp(x[s1], "c") && !strcmp(x[s2], "t"))
		Ns += 1;
	      else {
		if (!strcmp(x[s1], "t") && !strcmp(x[s2], "c"))
		  Ns += 1;
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L = *s;
      P = ((double) Ns/L);
      Q = ((double) (Nd - Ns)/L);
      a1 = 1 - 2 * P - Q;
      a2 = 1 - 2 * Q;
      if (*gamma) {
	b = -1 / *alpha;
	d[target] = *alpha * (pow(a1, b) + 0.5 * pow(a2, b) - 1.5) / 2;
      }
      else d[target] = -0.5 * log(a1 * sqrt(a2));
      if (*variance) {
	if (*gamma) {
	  b = -(1 / *alpha + 1);
	  c1 = pow(a1, b);
	  c2 = pow(a2, b);
	  c3 = (c1 + c2) / 2;
	} else {
	  c1 = 1 / a1;
	  c2 = 1 / a2;
	  c3 = (c1 + c2) / 2;
	}
	var[target] = (c1*c1*P + c3*c3*Q - pow(c1*P + c3*Q, 2))/L;
      }
    }
  }
}

void dist_dna_F81(char **x, int *n, int *s, double *d, double *BF,
		  int *pairdel, int *variance, double *var,
		  int *gamma, double *alpha)
{
  int i, j, k, Nd, L, target, s1, s2;
  double E, p;

  E = 1 - BF[0]*BF[0] - BF[1]*BF[1] - BF[2]*BF[2] - BF[3]*BF[3];

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) Nd += 1;
      }
      if (!*pairdel) L = *s;
      p = ((double) Nd/L);
      if (*gamma) d[target] = E * *alpha * (pow(1 - p/E, -1/ *alpha) - 1);
      else d[target] = -E*log(1 - p/E);
      if (*variance) {
	if (*gamma) var[target] = p*(1 - p)/(pow(1 - p/E, -2/(*alpha + 1)) * L);
	else var[target] = p*(1 - p)/(pow(1 - p/E, 2)*L);
      }
    }
  }
}

void dist_dna_K81(char **x, int *n, int *s, double *d, int *pairdel,
		   int *variance, double *var)
{
  int i, j, k, Nd, Nv1, Nv2, L, s1, s2, target;
  double P, Q, R, a1, a2, a3, a, b, c;

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Nv1 = 0;
      Nv2 = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a")) {
	    if (!strcmp(x[s2], "t")) Nv1 += 1;
	    if (!strcmp(x[s2], "c")) Nv2 += 1;
	  } else {
	    if (!strcmp(x[s1], "g")) {
	      if (!strcmp(x[s2], "c")) Nv1 += 1;
	      if (!strcmp(x[s2], "t")) Nv2 += 1;
	    } else {
	      if (!strcmp(x[s1], "c")) {
		if (!strcmp(x[j * *s + k], "g")) Nv1 += 1;
		if (!strcmp(x[j * *s + k], "a")) Nv2 += 1;
	      } else {
		if (!strcmp(x[s1], "t")) {
		  if (!strcmp(x[s2], "a")) Nv1 += 1;
		  if (!strcmp(x[s2], "g")) Nv2 += 1;
		}
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L = *s;
      P = ((double) (Nd - Nv1 - Nv2)/L);
      Q = ((double) Nv1/L);
      R = ((double) Nv2/L);
      a1 = 1 - 2*P - 2*Q;
      a2 = 1 - 2*P - 2*R;
      a3 = 1 - 2*Q - 2*R;
      d[target] = -0.25 * log(a1 * a2 * a3);
      if (*variance) {
        a = (1/a1 + 1/a2)/2;
	b = (1/a1 + 1/a3)/2;
	c = (1/a2 + 1/a3)/2;
        var[target] = (a*a*P + b*b*Q + c*c*R - pow(a*P + b*Q + c*R, 2))/2;
      }
    }
  }
}

void dist_dna_F84(char **x, int *n, int *s, double *d, double *BF,
		  int *pairdel, int *variance, double *var)
{
  int i, j, k, Nd, Ns, L, target, s1, s2;
  double P, Q, A, B, C, a, b, t1, t2, t3;

  A = (BF[0] * BF[2])/(BF[0] + BF[2]) + (BF[1] * BF[3])/(BF[1] + BF[3]);
  B = BF[0] * BF[2] + BF[1] * BF[3];
  C = (BF[0] + BF[2]) * (BF[1] + BF[3]);

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Ns = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a") && !strcmp(x[s2], "g"))
	    Ns += 1;
	  else {
	    if (!strcmp(x[s1], "g") && !strcmp(x[s2], "a"))
	      Ns += 1;
	    else {
	      if (!strcmp(x[s1], "c") && !strcmp(x[s2], "t"))
		Ns += 1;
	      else {
		if (!strcmp(x[s1], "t") && !strcmp(x[s2], "c"))
		  Ns += 1;
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L = *s;
      P = ((double) Ns/L);
      Q = ((double) (Nd - Ns)/L);
      d[target] = -2*A*log(1 - (P/(2*A) - (A - B)*Q/(2*A*C))) + 2*(A - B - C)*log(1 - Q/(2*C));
      if (*variance) {
        t1 = A*C;
	t2 = C*P/2;
	t3 = (A - B)*Q/2;
        a = t1/(t1 - t2 - t3);
	b = A*(A - B)/(t1 - t2 - t3) - (A - B - C)/(C - Q/2);
	var[target] = (a*a*P + b*b*Q - pow(a*P + b*Q, 2))/2;
      }
    }
  }
}

void dist_dna_T92(char **x, int *n, int *s, double *d, double *BF,
                  int *pairdel, int *variance, double *var)
{
  int i, j, k, Nd, Ns, L, target, s1, s2;
  double P, Q, wg, a1, a2, c1, c2, c3;

  wg = 2 * (BF[1] + BF[2]) * (1 - (BF[1] + BF[2]));

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Ns = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a") && !strcmp(x[s2], "g"))
	    Ns += 1;
	  else {
	    if (!strcmp(x[s1], "g") && !strcmp(x[s2], "a"))
	      Ns += 1;
	    else {
	      if (!strcmp(x[s1], "c") && !strcmp(x[s2], "t"))
		Ns += 1;
	      else {
		if (!strcmp(x[s1], "t") && !strcmp(x[s2], "c"))
		  Ns += 1;
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L = *s;
      P = ((double) Ns/L);
      Q = ((double) (Nd - Ns)/L);
      a1 = 1 - P/wg - Q;
      a2 = 1 - 2*Q;
      d[target] = -wg*log(a1) - 0.5*(1 - wg)*log(a2);
      if (*variance) {
        c1 = 1/a1;
        c2 = 1/a2;
        c3 = wg*(c1 - c2) + c2;
        var[target] = (c1*c1*P + c3*c3*Q - pow(c1*P + c3*Q, 2))/L;
      }
    }
  }
}

void dist_dna_TN93(char **x, int *n, int *s, double *d, double *BF,
		  int *pairdel, int *variance, double *var,
		  int *gamma, double *alpha)
{
  int i, j, k, Nd, Ns1, Ns2, L, target, s1, s2;
  double P1, P2, Q, A, B, C, gR, gY, k1, k2, k3, k4, w1, w2, w3, c1, c2, c3, c4, b;

  gR = BF[0] + BF[2];
  gY = BF[1] + BF[3];
  k1 = 2 * BF[0] * BF[2] / gR;
  k2 = 2 * BF[1] * BF[3] / gY;
  k3 = 2 * (gR * gY - BF[0]*BF[2]*gY/gR - BF[1]*BF[3]*gR/gY);

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Ns1 = 0;
      Ns2 = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
        s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a") && !strcmp(x[s2], "g"))
	    Ns1 += 1;
	  else {
	    if (!strcmp(x[s1], "g") && !strcmp(x[s2], "a"))
	      Ns1 += 1;
	    else {
	      if (!strcmp(x[s1], "c") && !strcmp(x[s2], "t"))
		Ns2 += 1;
	      else {
		if (!strcmp(x[s1], "t") && !strcmp(x[s2], "c"))
		  Ns2 += 1;
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L = *s;
      P1 = ((double) Ns1/L);
      P2 = ((double) Ns2/L);
      Q = ((double) (Nd - Ns1 - Ns2)/L);
      w1 = 1 - P1/k1 - Q/(2*gR);
      w2 = 1 - P2/k2 - Q/(2*gY);
      w3 = 1 - Q/(2*gR*gY);
      if (*gamma) {
        k4 = 2*(BF[0]*BF[2] + BF[1]*BF[3] + gR*gY);
        b = -1 / *alpha;
        c1 = pow(w1, b);
        c2 = pow(w2, b);
        c3 = pow(w3, b);
        c4 = k1*c1/(2*gR) + k2*c2/(2*gY) + k3*c3/(2*gR*gY);
        d[target] = *alpha * (k1*pow(w1, b) + k2*pow(w2, b) + k3*pow(w3, b) - k4);
      } else {
        k4 = 2*((BF[0]*BF[0] + BF[2]*BF[2])/(2*gR*gR) + (BF[2]*BF[2] + BF[3]*BF[3])/(2*gY*gY));
        c1 = 1/w1;
        c2 = 1/w2;
        c3 = 1/w3;
        c4 = k1 * c1/(2 * gR) + k2 * c2/(2 * gY) + k4 * c3;
        d[target] = -k1*log(w1) - k2*log(w2) - k3*log(w3);
      }
      if (*variance) {
        var[target] = (c1*c1*P1 + c2*c2*P2 + c4*c4*Q - pow(c1*P1 + c2*P2 + c4*Q, 2))/L;
      }
    }
  }
}

void dist_dna_GG95(char **x, int *n, int *s, double *d,
		   int *pairdel, int *variance, double *var)
{
  int i, j, k, s1, s2, target, *GC, gccount, Nd, Ns, *L, length, tl, npair;
  double *theta, gcprop, *P, pp, *Q, qq, *tstvr, svr, A,
    sum, ma /* mean alpha */, K1, K2;

  GC = &gccount;
  theta = &gcprop;
  L = &length;
  P = &pp;
  Q = &qq;
  tstvr = &svr;

  npair = *n * (*n - 1) / 2;

  GC = (int*)R_alloc(*n, sizeof(int));
  theta = (double*)R_alloc(*n, sizeof(double));
  L = (int*)R_alloc(npair, sizeof(int));
  P = (double*)R_alloc(npair, sizeof(double));
  Q = (double*)R_alloc(npair, sizeof(double));
  tstvr = (double*)R_alloc(npair, sizeof(double));

  /* get the proportion of GC, theta, in each sequence */
  for (i = 0; i < *n; i++) {
    tl = 0;
    GC[i] = 0;
    for (k = 0; k < *s; k++) {
      s1 = i * *s + k;
      if (*pairdel) {
	if (strcmp(x[s1], "n")) tl += 1;
	else continue;
      }
      if (!strcmp(x[s1], "c") || !strcmp(x[s1], "g"))
	GC[i] += 1;
    }
    if (!*pairdel) tl = *s;
    theta[i] = ((double) GC[i]/tl);
  }

  /* get the proportions of transitions and transversions,
     and the estimates of their ratio for each pair; we
     also get the sample size for each pair in L  */
  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Ns = 0;
      L[target] = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L[target] += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a") && !strcmp(x[s2], "g"))
	    Ns += 1;
	  else {
	    if (!strcmp(x[s1], "g") && !strcmp(x[s2], "a"))
	      Ns += 1;
	    else {
	      if (!strcmp(x[s1], "c") && !strcmp(x[s2], "t"))
		Ns += 1;
	      else {
		if (!strcmp(x[s1], "t") && !strcmp(x[s2], "c"))
		  Ns += 1;
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L[target] = *s;
      P[target] = ((double) Ns/L[target]);
      Q[target] = ((double) (Nd - Ns)/L[target]);
      A = log(1 - 2*Q[target]);
      tstvr[target] = 2*(log(1 - 2*P[target] - Q[target]) - 0.5*A)/A;
    }
  }

  /* compute the mean alpha (ma) = mean Ts/Tv */
  sum = 0;
  tl = 0;
  for (i = 0; i < npair; i++)
    /* some values of tstvr are -inf if there is no
       transversions observed: we exclude them */
    if (R_FINITE(tstvr[i])) {
      sum += tstvr[i];
      tl += 1;
    }
  ma = sum/tl;

  /* compute the distance for each pair */
  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      A = 1 - 2*Q[target];
      K1 = 1 + ma*(theta[i]*(1 - theta[i]) + theta[j]*(1 - theta[j]));
      K2 = ma*pow(theta[i] - theta[j], 2)/(ma + 1);
      d[target] = -0.5*K1*log(A) + K2*(1 - pow(A, 0.25*(ma + 1)));
      if (*variance)
	var[target] = pow(K1 + K2*0.5*(ma + 1)*pow(A, 0.25*(ma + 1)), 2)*Q[target]*(1 - Q[target])/(A*A*L[target]);
    }
  }
}

void dist_dna(char **x, int *n, int *s, int *model, double *d,
	      double *BF, int *pairdel, int *variance, double *var,
	      int *gamma, double *alpha)
{
  switch (*model) {
  case 1 : dist_dna_JC69(x, n, s, d, pairdel, variance,
			 var, gamma, alpha); break;
  case 2 : dist_dna_K80(x, n, s, d, pairdel, variance,
			var, gamma, alpha); break;
  case 3 : dist_dna_F81(x, n, s, d, BF, pairdel, variance,
			var, gamma, alpha); break;
  case 4 : dist_dna_K81(x, n, s, d, pairdel, variance,
			var); break;
  case 5 : dist_dna_F84(x, n, s, d, BF, pairdel, variance,
			var); break;
  case 6 : dist_dna_T92(x, n, s, d, BF, pairdel, variance,
			var); break;
  case 7 : dist_dna_TN93(x, n, s, d, BF, pairdel, variance,
			 var, gamma, alpha); break;
  case 8 : dist_dna_GG95(x, n, s, d, pairdel,
			 variance, var); break;
  }
}
