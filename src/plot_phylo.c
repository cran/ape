/* Copyright 2004 Emmanuel Paradis <paradis@isem.univ-montp2.fr> */

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

void node_depth_edgelength(int *ntip, int *nnode, int *edge1, int *edge2, int *nms, double *edge_length, double *xx)
{
  int i, j, k;

  for (i = 1; i < *ntip + *nnode; i++) {
    j = 0;
    while (edge2[j] != nms[i]) j++;
    if (edge1[j] < 0) k = -edge1[j] - 1; 
    else k = nnode + edge1[j] - 1;
    xx[i] = xx[k] + edge_length[j];
  }
}
