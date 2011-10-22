/* ape.c    2011-10-11 */

/* Copyright 2011 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "ape.h"

int give_index(int i, int j, int n)
{
	if (i > j) return(DINDEX(j, i));
	else return(DINDEX(i, j));
}
