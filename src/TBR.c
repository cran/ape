/* TBR.c    2014-03-04 */

/* Copyright 2009 Richard Desper */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "me.h"

/*functions from me_balanced.c*/
void makeBMEAveragesTable(tree *T, double **D, double **A);
void assignBMEWeights(tree *T, double **A);

/*from me.c*/
edge *siblingEdge(edge *e);
double **initDoubleMatrix(int d);
void freeMatrix(double **D, int size);
edge *depthFirstTraverse(tree *T, edge *e);
edge *findBottomLeft(edge *e);

/*from bnni.c*/
void weighTree(tree *T);

void freeTree(tree *T);

/*from SPR.c*/
void zero3DMatrix(double ***X, int h, int l, int w);

void assignTBRDownWeightsUp(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights,
			    double *bestWeight, edge **eBestSplit, edge **eBestTop, edge **eBestBottom);
void assignTBRDownWeightsSkew(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights,
			      double *bestWeight, edge **eBestSplit, edge **eBestTop, edge **eBestBototm);
void assignTBRDownWeightsDown(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights,
			      double *bestWeight, edge **eBestSplit, edge **eBestTop, edge **eBestBottom);

void assignTBRUpWeights(edge *ebottom, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A,
						double **dXTop, double ***swapWeights, edge *etop, double *bestWeight,
						edge **bestSplit, edge **bestTop, edge **bestBottom)
						/*function assigns the value for etop if the tree below vtest is moved to be below etop*/
						/*and SPR for tree bottom tree splits ebottom*/
						/*recursive function searching over values of ebottom*/
						/*minor variant of SPR.c's assignUpWeights
						difference is the index assignment in the array swapWeights, which has a different meaning
						for the TBR routines*/
/*also using dXTop to assign value of average distance to tree above vtest*/
{ /*SPR performed on tree above vtest...*/
  edge *sib, *left, *right;
  /*B is above vtest, A is other tree below vtest unioned with trees in path to vtest*/
  /*sib is tree C being passed by B*/
  /*D is tree below etest*/
  double D_AB, D_CD, D_AC, D_BD;
  sib = siblingEdge(ebottom);
  left = ebottom->head->leftEdge;
  right = ebottom->head->rightEdge;
  if (NULL == etop)
    {
      if (NULL == back) /*first recursive call*/
	{
	  if (NULL == left)
	    return; /*no subtree below for SPR*/
	  else       /*NULL == back and NULL == etop*/
	    {
	      assignTBRUpWeights(left,vtest,va,ebottom,va,A[va->index][vtest->index],0.5,A,dXTop,swapWeights,NULL,bestWeight,bestSplit,bestTop,bestBottom);
	      assignTBRUpWeights(right,vtest,va,ebottom,va,A[va->index][vtest->index],0.5,A,dXTop,swapWeights,NULL,bestWeight,bestSplit,bestTop,bestBottom);
	    }
	}
      else /*NULL != back*/
	{
	  D_BD = A[ebottom->head->index][vtest->index];
	  /*B is tree above vtest, D is tree below ebottom*/
	  D_CD = A[sib->head->index][ebottom->head->index]; /*C is tree below sib*/
	  D_AC = A[back->head->index][sib->head->index] +
	    coeff*(A[va->index][sib->head->index] - A[vtest->index][sib->head->index]);
	  /*va is root of subtree skew back at vtest*/
	  /*A is union of tree below va and all subtrees already passed in path from vtest to ebottom*/
	  D_AB = 0.5*(oldD_AB + A[vtest->index][cprev->index]);
	  swapWeights[vtest->index][ebottom->head->index][ebottom->head->index] = swapWeights[vtest->index][back->head->index][back->head->index] + (D_AC + D_BD - D_AB - D_CD);
	  if (swapWeights[vtest->index][ebottom->head->index][ebottom->head->index] < *bestWeight)
	    {
	      *bestSplit = vtest->parentEdge;
	      *bestTop = NULL;
	      *bestBottom = ebottom;
	      *bestWeight = swapWeights[vtest->index][ebottom->head->index][ebottom->head->index];
	    }
	  if (NULL != left)
	    {
	      assignTBRUpWeights(left,vtest,va,ebottom,sib->head,D_AB,0.5*coeff,A,dXTop,swapWeights,NULL,bestWeight,bestSplit,bestTop,bestBottom);
	      assignTBRUpWeights(right,vtest,va,ebottom,sib->head,D_AB,0.5*coeff,A,dXTop,swapWeights,NULL,bestWeight,bestSplit,bestTop,bestBottom);
	    }
	}
    }
  else /*NULL != etop*/
    {
      if (NULL == back) /*first recursive call*/
	{
	  if (swapWeights[vtest->index][etop->head->index][etop->head->index]< *bestWeight)
	    /*this represents value of SPR from esplit to etop, with no SPR in bottom tree*/
	    {
	      *bestSplit = vtest->parentEdge;
	      *bestTop = etop;
	      *bestBottom = NULL;
	      *bestWeight = swapWeights[vtest->index][etop->head->index][etop->head->index];
	    }
	  if (NULL == left)
	    return; /*no subtree below for SPR*/
	  else if (NULL != etop)/*start the process of assigning weights recursively*/
	    {
	      assignTBRUpWeights(left,vtest,va,ebottom,va,dXTop[va->index][etop->head->index],0.5,A,dXTop,swapWeights,etop,bestWeight,bestSplit,bestTop,bestBottom);
	      assignTBRUpWeights(right,vtest,va,ebottom,va,dXTop[va->index][etop->head->index],0.5,A,dXTop,swapWeights,etop,bestWeight,bestSplit,bestTop,bestBottom);
	    }
	} /*NULL == back*/
      /*in following bit, any average distance of form A[vtest->index][x->index] is
	replaced by dXTop[x->index][etop->head->index]*/
      else /*second or later recursive call, NULL != etop*/
	{
	  D_BD = dXTop[ebottom->head->index][etop->head->index]; /*B is tree above vtest - it is in configuration
								   indexed by etop*/
	  /*D is tree below ebottom*/
	  D_CD = A[sib->head->index][ebottom->head->index]; /*C is tree below sib*/
	  D_AC = A[back->head->index][sib->head->index] +
	    coeff*(A[va->index][sib->head->index] - A[sib->head->index][vtest->index]);
	  /*it is correct to use A[][] here because the bad average distances involving B from the first term will
	    be cancelled by the bad average distances involving B in the subtracted part*/
	  /*va is root of subtree skew back at vtest*/
	  /*A is union of tree below va and all subtrees already passed in path from vtest to ebottom*/
	  D_AB = 0.5*(oldD_AB + dXTop[cprev->index][etop->head->index]);
	  swapWeights[vtest->index][ebottom->head->index][etop->head->index] = swapWeights[vtest->index][back->head->index][etop->head->index]  + (D_AC + D_BD - D_AB - D_CD);
	  if (swapWeights[vtest->index][ebottom->head->index][etop->head->index] + swapWeights[vtest->index][etop->head->index][etop->head->index]< *bestWeight)
	    /*first term is contribution of second SPR, second term is contribution of first SPR*/
	    {
	      *bestSplit = vtest->parentEdge;
	      *bestTop = etop;
	      *bestBottom = ebottom;
	      *bestWeight = swapWeights[vtest->index][ebottom->head->index][etop->head->index] + swapWeights[vtest->index][etop->head->index][etop->head->index];
	    }
	  if (NULL != left)
	    {
	      assignTBRUpWeights(left,vtest, va, ebottom, sib->head, D_AB, 0.5*coeff, A, dXTop, swapWeights,etop,bestWeight,bestSplit,bestTop,bestBottom);
	      assignTBRUpWeights(right,vtest, va, ebottom, sib->head, D_AB, 0.5*coeff, A, dXTop, swapWeights,etop,bestWeight,bestSplit,bestTop,bestBottom);
	    }
	} /*else NULL != back, etop*/
    }
}

/*recall NNI formula: change in tree length from AB|CD split to AC|BD split is
proportional to D_AC + D_BD - D_AB - D_CD*/
/*in our case B is the tree being moved (below vtest), A is the tree backwards below back, but
with the vtest subtree removed, C is the sibling tree of back and D is the tree above vtest*/
/*use va to denote the root of the sibling tree to B in the original tree*/
/*please excuse the multiple uses of the same letters: A,D, etc.*/
void assignTBRDownWeightsUp(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights,
			    double *bestWeight, edge **bestSplitEdge, edge **bestTop, edge **bestBottom)
{
	edge *par, *sib, *skew;
	double D_AC, D_BD, D_AB, D_CD;
	par = etest->tail->parentEdge;
	skew = siblingEdge(etest);
	if (NULL == back) /*first recursive call*/
	  {
	  if (NULL == par)
	    return;
	  else /*start the process of assigning weights recursively*/
	    {
	      assignTBRDownWeightsUp(par,vtest,va,etest,va,A[va->index][vtest->index],0.5,A,swapWeights,bestWeight,bestSplitEdge,bestTop,bestBottom);
	      assignTBRDownWeightsSkew(skew,vtest,va,etest,va,A[va->index][vtest->index],0.5,A,swapWeights,bestWeight,bestSplitEdge,bestTop,bestBottom);
	    }
	}
	else /*second or later recursive call*/
	  {
	    sib = siblingEdge(back);
	    D_BD = A[vtest->index][etest->head->index]; /*straightforward*/
	    D_CD = A[sib->head->index][etest->head->index]; /*this one too*/
	    D_AC = A[sib->head->index][back->head->index] + coeff*(A[sib->head->index][va->index]
								   - A[sib->head->index][vtest->index]);
	    D_AB = 0.5*(oldD_AB + A[vtest->index][cprev->index]);
	    swapWeights[vtest->index][etest->head->index][etest->head->index] = swapWeights[vtest->index][back->head->index][back->head->index] + (D_AC + D_BD - D_AB - D_CD);
	    /*using diagonal to store values for SPR swaps above the split edge*/
	    /*this is value of swapping tree below vtest to break etest*/
	    if (swapWeights[vtest->index][etest->head->index][etest->head->index] < *bestWeight)
	      {
		*bestWeight = swapWeights[vtest->index][etest->head->index][etest->head->index];
		*bestSplitEdge = vtest->parentEdge;
		*bestTop = etest;
		*bestBottom = NULL;
	      }
	    if (NULL != par)
	      {
		assignTBRDownWeightsUp(par,vtest,va,etest,sib->head,D_AB,0.5*coeff,A,swapWeights,bestWeight,bestSplitEdge,bestTop,bestBottom);
		assignTBRDownWeightsSkew(skew,vtest,va,etest,sib->head,D_AB,0.5*coeff,A,swapWeights,bestWeight,bestSplitEdge,bestTop,bestBottom);
	      }
	  }
}


void assignTBRDownWeightsSkew(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights,
			      double *bestWeight, edge **bestSplitEdge, edge **bestTop, edge **bestBottom)
{
  /*same general idea as assignDownWeights, except needing to keep track of things a bit differently*/
  edge *par, *left, *right;
  /*par here = sib before
    left, right here = par, skew before*/
  double D_AB, D_CD, D_AC, D_BD;
  /*B is subtree being moved - below vtest
    A is subtree remaining fixed - below va, unioned with all trees already passed by B*/
  /*C is subtree being passed by B, in this case above par
    D is subtree below etest, fixed on other side*/
  par = etest->tail->parentEdge;
  left = etest->head->leftEdge;
  right = etest->head->rightEdge;
  if (NULL == back)
    {
      if (NULL == left)
	return;
      else
	{
	  assignTBRDownWeightsDown(left,vtest,va,etest,etest->tail,A[vtest->index][etest->tail->index],0.5,A,swapWeights,bestWeight,bestSplitEdge,bestTop,bestBottom);
	  assignTBRDownWeightsDown(right,vtest,va,etest,etest->tail,A[vtest->index][etest->tail->index],0.5,A,swapWeights,bestWeight,bestSplitEdge,bestTop,bestBottom);
	}
    }
  else
    {
      D_BD = A[vtest->index][etest->head->index];
      D_CD = A[par->head->index][etest->head->index];
      D_AC = A[back->head->index][par->head->index] + coeff*(A[va->index][par->head->index] - A[vtest->index][par->head->index]);
      D_AB = 0.5*(oldD_AB + A[vtest->index][cprev->index]);
      swapWeights[vtest->index][etest->head->index][etest->head->index] = swapWeights[vtest->index][back->head->index][back->head->index] + (D_AC + D_BD - D_AB - D_CD);
      if (swapWeights[vtest->index][etest->head->index][etest->head->index] < *bestWeight)
	{
	  *bestWeight = swapWeights[vtest->index][etest->head->index][etest->head->index];
	  *bestSplitEdge = vtest->parentEdge;
	  *bestTop = etest;
	  *bestBottom = NULL;
	}
      if (NULL != left)
	{
	  assignTBRDownWeightsDown(left,vtest, va, etest, etest->tail, D_AB, 0.5*coeff, A, swapWeights,bestWeight,bestSplitEdge,bestTop,bestBottom);
	  assignTBRDownWeightsDown(right,vtest, va, etest, etest->tail, D_AB, 0.5*coeff, A, swapWeights,bestWeight,bestSplitEdge,bestTop,bestBottom);
	}
    }
}

void assignTBRDownWeightsDown(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights,
			      double *bestWeight, edge **bestSplitEdge, edge **bestTop, edge **bestBottom)
{
  /*again the same general idea*/
  edge *sib, *left, *right;
  /*sib here = par in assignDownWeightsSkew
    rest is the same as assignDownWeightsSkew*/
  double D_AB, D_CD, D_AC, D_BD;
  /*B is below vtest, A is below va unioned with all trees already passed by B*/
  /*C is subtree being passed - below sib*/
  /*D is tree below etest*/
  sib = siblingEdge(etest);
  left = etest->head->leftEdge;
  right = etest->head->rightEdge;
  D_BD = A[vtest->index][etest->head->index];
  D_CD = A[sib->head->index][etest->head->index];
  D_AC = A[sib->head->index][back->head->index] + coeff*(A[sib->head->index][va->index] - A[sib->head->index][vtest->index]);
  D_AB = 0.5*(oldD_AB + A[vtest->index][cprev->index]);
  swapWeights[vtest->index][etest->head->index][etest->head->index] = swapWeights[vtest->index][back->head->index][back->head->index] + ( D_AC + D_BD - D_AB - D_CD);
  if (swapWeights[vtest->index][etest->head->index][etest->head->index] < *bestWeight)
    {
      *bestWeight = swapWeights[vtest->index][etest->head->index][etest->head->index];
      *bestSplitEdge = vtest->parentEdge;
      *bestTop = etest;
      *bestBottom = NULL;
    }
  if (NULL != left)
    {
      assignTBRDownWeightsDown(left,vtest, va, etest, sib->head, D_AB, 0.5*coeff, A, swapWeights,bestWeight,bestSplitEdge,bestTop,bestBottom);
      assignTBRDownWeightsDown(right,vtest, va, etest, sib->head, D_AB, 0.5*coeff, A, swapWeights,bestWeight,bestSplitEdge,bestTop,bestBottom);
    }
}

/*general idea is to have a double loop for a given edge, testing all SPRs for the subtrees above and below a given edge.
  Then that function loops over all the edges of a tree*/

void TBRswitch(tree *T, edge *e1, edge *e2, edge *e3);

/*vbottom is node below esplit for average calculations in matrix dXTop, A is matrix of average
  distances from original tree, dXout is average distance from vbottom to tree rooted at far edge
  of eback, if SPR breaking eback, UpOrDown indicates whether etop is in path above split edge
  (Up) or not (Down)*/
void calcTBRTopBottomAverage(node *vbottom, double **A, double **dXTop, double dXOut,
			     edge *esplit, edge *etop, edge *eback, int UpOrDown)
{
  edge *enew1, *enew2, *emove;
  double newdXOut;
  if (esplit == eback) /*top level call - means trivial SPR*/
    dXTop[vbottom->index][etop->head->index] = A[vbottom->index][esplit->head->index];
  else
    dXTop[vbottom->index][etop->head->index] = dXTop[vbottom->index][eback->head->index] +
      0.25*(A[vbottom->index][etop->head->index] - dXOut);
  /*by moving etop past the vbottom tree, everything in the etop tree is closer by coefficient of
    0.25, while everything in the old back tree is further by a coefficient of 0.25*/
  /*everything in the tree that is being moved (emove) keeps the same relative weight in the average
    distance calculation*/
  if (UP == UpOrDown)
    {
      enew1 = etop->tail->parentEdge;
      if (NULL != enew1) /*need to do recursive calls*/
	{
	  enew2 = siblingEdge(etop);
	  emove = siblingEdge(eback);	 /*emove is third edge meeting at vertex with eback, etest*/
	  if (esplit == eback)
	    newdXOut = A[vbottom->index][emove->head->index];
	  else
	    newdXOut = 0.5*(dXOut + A[vbottom->index][emove->head->index]);
	  calcTBRTopBottomAverage(vbottom,A,dXTop,newdXOut,esplit, enew1,etop,UP); /*old etop is new value for back*/
	  calcTBRTopBottomAverage(vbottom,A,dXTop,newdXOut,esplit, enew2,etop,DOWN);
	}
    }
  else /*moving down*/
    {
      enew1 = etop->head->leftEdge;
      if (NULL != enew1)
	{
	  enew2 = etop->head->rightEdge;
	  if (eback == siblingEdge(etop))
	    emove = etop->tail->parentEdge;
	  else
	    emove = siblingEdge(etop);
	  if (esplit == eback)
	    newdXOut = A[vbottom->index][emove->head->index];
	  else
	    newdXOut = 0.5*(dXOut + A[vbottom->index][emove->head->index]);
	  calcTBRTopBottomAverage(vbottom,A,dXTop,newdXOut,esplit,enew1,etop,DOWN);
	  calcTBRTopBottomAverage(vbottom,A,dXTop,newdXOut,esplit,enew2,etop,DOWN);
	}
    }
}

void calcTBRaverages(tree *T, edge *esplit, double **A, double **dXTop)
{
  edge *ebottom, *par, *sib;
  for (ebottom = findBottomLeft(esplit); ebottom != esplit; ebottom = depthFirstTraverse(T,ebottom))
    {
      par = esplit->tail->parentEdge;
      sib = siblingEdge(esplit);
      calcTBRTopBottomAverage(ebottom->head,A, dXTop, 0.0, esplit, par,esplit,UP);
      calcTBRTopBottomAverage(ebottom->head,A, dXTop, 0.0, esplit, sib,esplit,DOWN);
    }
}

void TBR(tree *T, double **D, double **A)
{
  int i;
  edge *esplit, *etop, *eBestTop, *eBestBottom, *eBestSplit;
  edge *eout, *block;
  edge *left, *right, *par, *sib;
  double **dXTop; /*dXTop[i][j] is average distance from subtree rooted at i to tree above split edge, if
		    SPR above split edge cuts edge whose head has index j*/
  double bestWeight;
  double ***TBRWeights;
  dXTop = initDoubleMatrix(T->size);
  weighTree(T);
  TBRWeights = (double ***)calloc(T->size,sizeof(double **));
  for(i=0;i<T->size;i++)
    TBRWeights[i] = initDoubleMatrix(T->size);
  do
    {
      zero3DMatrix(TBRWeights,T->size,T->size,T->size);
      bestWeight = 0.0;
      eBestSplit = eBestTop = eBestBottom = NULL;
      for(esplit=depthFirstTraverse(T,NULL);NULL!=esplit;esplit=depthFirstTraverse(T,esplit))
	{
	  par = esplit->tail->parentEdge;
	  if (NULL != par)
	    {
	      sib = siblingEdge(esplit);
	      /*next two lines calculate the possible improvements for any SPR above esplit*/
	      assignTBRDownWeightsUp(par,esplit->head,sib->head,NULL,NULL,0.0,1.0,A,TBRWeights,&bestWeight,&eBestSplit,&eBestTop,&eBestBottom);
	      assignTBRDownWeightsSkew(sib,esplit->head,sib->tail,NULL,NULL,0.0,1.0,A,TBRWeights,&bestWeight,&eBestSplit,&eBestTop,&eBestBottom);
	      calcTBRaverages(T,esplit,A,dXTop); /*calculates the average distance from any subtree
						   below esplit to the entire subtree above esplit,
						   after any possible SPR above*/
	      /*for etop above esplit, we loop using information from above to calculate values
		for all possible SPRs below esplit*/
	    }

	  right = esplit->head->rightEdge;
	  if (NULL != right)
	    {
	      left = esplit->head->leftEdge;
	      /*test case: etop = null means only do bottom SPR*/
	      assignTBRUpWeights(left,esplit->head,right->head,NULL,NULL,0.0,1.0,A,dXTop,TBRWeights,NULL,&bestWeight,&eBestSplit,&eBestTop,&eBestBottom);
	      assignTBRUpWeights(right,esplit->head,left->head,NULL,NULL,0.0,1.0,A,dXTop,TBRWeights,NULL,&bestWeight,&eBestSplit,&eBestTop,&eBestBottom);

	      block = esplit;
	      while (NULL != block)
		{
		  if (block != esplit)
		    {
		      etop = block;
		      assignTBRUpWeights(left,esplit->head,right->head,NULL,NULL,0.0,1.0,A,dXTop,TBRWeights,etop,&bestWeight,&eBestSplit,&eBestTop,&eBestBottom);
		      assignTBRUpWeights(right,esplit->head,left->head,NULL,NULL,0.0,1.0,A,dXTop,TBRWeights,etop,&bestWeight,&eBestSplit,&eBestTop,&eBestBottom);
		    }
		  eout = siblingEdge(block);
		  if (NULL != eout)
		    {
		      etop = findBottomLeft(eout);
		      while (etop->tail != eout->tail)
			{
			  /*for ebottom below esplit*/

			  assignTBRUpWeights(left,esplit->head,right->head,NULL,NULL,0.0,1.0,A,dXTop,TBRWeights,etop,&bestWeight,&eBestSplit,&eBestTop,&eBestBottom);
			  assignTBRUpWeights(right,esplit->head,left->head,NULL,NULL,0.0,1.0,A,dXTop,TBRWeights,etop,&bestWeight,&eBestSplit,&eBestTop,&eBestBottom);
			  etop = depthFirstTraverse(T,etop);
			}

		      /*etop == eout*/

		      assignTBRUpWeights(left,esplit->head,right->head,NULL,NULL,0.0,1.0,A,dXTop,TBRWeights,etop,&bestWeight,&eBestSplit,&eBestTop,&eBestBottom);
		      assignTBRUpWeights(right,esplit->head,left->head,NULL,NULL,0.0,1.0,A,dXTop,TBRWeights,etop,&bestWeight,&eBestSplit,&eBestTop,&eBestBottom);
		    }
		  block = block->tail->parentEdge;
		}
	    } /*if NULL != right*/
	} /*for esplit*/
      /*find bestWeight, best split edge, etc.*/
      if (bestWeight < -EPSILON)
	{
//	  if (verbose)
//	    {
//	      printf("TBR #%d: Splitting edge %s: top edge %s, bottom edge %s\n",*count+1,
//		     eBestSplit->label, eBestTop->label,eBestBottom->label);
//	      printf("Old tree weight is %lf, new tree weight should be %lf\n",T->weight, T->weight + 0.25*bestWeight);
//	    }
	  TBRswitch(T,eBestSplit,eBestTop,eBestBottom);
	  makeBMEAveragesTable(T,D,A);
	  assignBMEWeights(T,A);
	  weighTree(T);
//	  if (verbose)
//	    printf("TBR: new tree weight is %lf\n\n",T->weight);<
//	  (*count)++;
	}
      else
	bestWeight = 1.0;
    } while (bestWeight < -EPSILON);
  for(i=0;i<T->size;i++)
    freeMatrix(TBRWeights[i],T->size);
  freeMatrix(dXTop,T->size);
  free(TBRWeights); /* added by EP 2014-03-04 */
}

void SPRTopShift(tree *T, node *v, edge *e, int UpOrDown);

void TBRswitch(tree *T, edge *es, edge *et, edge *eb)
{
	if (NULL != et)
		SPRTopShift(T,es->head,et,DOWN); /*DOWN because tree being moved is below split edge*/
	if (NULL != eb)
		SPRTopShift(T,es->head,eb,UP);   /*UP because tree being moved is above split edge*/
}
