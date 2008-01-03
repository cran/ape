/* 
 *  nprsfunc.c
 *
 * (c) 2003  Gangolf Jobb and Korbinian Strimmer
 *
 *  Functions for nonparametric rate smoothing (NPRS)
 *  (see MJ Sanderson. 1997.  MBE 14:1218-1231)
 *
 *  This code may be distributed under the GNU GPL
 */


#include <stdlib.h>
#include <string.h>
#include <math.h>

/*============================================================================*/

#define EPS 1e-6
#define LARGENUM 1e99

/* original scale */
#define MINPARAM EPS
#define MAXPARAM 1-EPS
#define STARTPARAM 0.5


/* log scale */
#define MINLPARAM log(MINPARAM)
#define MAXLPARAM log(MAXPARAM)
#define STARTLPARAM log(STARTPARAM)

#define ARRLEN 2048
#define LABLEN 64

#define index aperates_index

int tree_lowerNodes[ARRLEN];
int tree_upperNodes[ARRLEN];
double tree_edgeLengths[ARRLEN];
int tree_nedges;
char tree_tipLabels[ARRLEN][LABLEN];
int tree_ntips;

int nparams,index[ARRLEN],ancestor[ARRLEN];

/*============================================================================*/

void setTree(
 int *lowerNodes,
 int *upperNodes,
 double *edgeLengths,
 double *minEdgeLength,
 int *nedges,
 char **tipLabels,
 int *ntips,
 int *result
) {
 int i;

 tree_nedges=*nedges;

 for(i=0;i<tree_nedges;i++) {
  tree_lowerNodes[i]=lowerNodes[i];
  tree_upperNodes[i]=upperNodes[i];
  tree_edgeLengths[i]=(edgeLengths[i]<*minEdgeLength)?*minEdgeLength:edgeLengths[i];
 }

 tree_ntips=*ntips;

 for(i=0;i<tree_ntips;i++) {
  strcpy(tree_tipLabels[i],tipLabels[i]);
 }

 nparams=0;
 for(i=0;i<tree_nedges;i++) {
  if(tree_lowerNodes[i]<0&&tree_upperNodes[i]<0) {index[-tree_upperNodes[i]]=nparams++;}
  if(tree_upperNodes[i]<0) ancestor[-tree_upperNodes[i]]=tree_lowerNodes[i];   
 }

 *result=0;
}

/*============================================================================*/

void getNFreeParams(int *result) { *result=nparams; }

/*============================================================================*/

void getNEdges(int *result) { *result=tree_nedges; }

/*============================================================================*/

void getEdgeLengths(double *result) {
 int i;

 for(i=0;i<tree_nedges;i++) result[i]=tree_edgeLengths[i];

}

/*============================================================================*/

double age(double *params,int node) {
 double prod;

 if(node>=0) return(0.0);
 if(node==-1) return(1.0);

 prod=0.0;
 while(node!=-1) {prod+=params[index[-node]]; node=ancestor[-node];}

 return(exp(prod)); 
}


void getDurations(double *params,double *scale,double *result) {
 int i,low,upp;

 for(i=0;i<tree_nedges;i++) {
  low=tree_lowerNodes[i]; upp=tree_upperNodes[i];
  if(low<0&&upp<0)  result[i]=*scale*(age(params,low)-age(params,upp));
  else              result[i]=*scale*age(params,low);
 }

}

/*============================================================================*/

void getRates(double *params,double *scale,double *result) {
 int i,low,upp;

 for(i=0;i<tree_nedges;i++) {
  low=tree_lowerNodes[i]; upp=tree_upperNodes[i];
  if(low<0&&upp<0)  result[i]=tree_edgeLengths[i]/(*scale*(age(params,low)-age(params,upp)));
  else              result[i]=tree_edgeLengths[i]/(*scale*age(params,low));
 }

}

/*============================================================================*/


/* note on the choice of parameters:
 * 
 * in order to obtain a chronogram we optimize the
 * relative heights of each internal node, i.e. to
 * each internal node (minus the root) we assign a number
 * between 0 and 1 (MINPARAM - MAXPARAM).
 *
 * To obtain the actual height of a node the relative heights
 * of the nodes have to multiplied (in fact we work on the log-scale
 * so we simply sum).
 */


/* internal objective function - parameters are on log scale */
void nprsObjectiveFunction(
 double *params,
 int *expo,
 double *result
) {
 int p,k,j;
 double me,rk,rj,sum,scale=1.0,durations[ARRLEN],rates[ARRLEN];

 p=*expo;
                  
 getDurations(params,&scale,durations);

 for(k=0;k<tree_nedges;k++) rates[k]=tree_edgeLengths[k]/durations[k];

 me=0.; for(k=0;k<tree_nedges;k++) me+=rates[k]; me/=tree_nedges;

 sum=0.;

 for(k=0;k<tree_nedges;k++) {                      rk=rates[k];
  for(j=0;j<tree_nedges;j++) { if(j==k) continue;  rj=rates[j];
   if(tree_lowerNodes[j]==-1)                  sum+=pow(fabs(me-rj),p);  else
   if(tree_upperNodes[k]==tree_lowerNodes[j])  sum+=pow(fabs(rk-rj),p);
  }
 }

 *result=sum;

}


/* check parameter bounds on log scale */
int checkLogParams(double *params)
{
   int i;
   for(i=0; i<nparams; i++)
   {
     if(params[i] > MAXLPARAM || params[i] < MINLPARAM) return 0;  /* out of bounds */
   }
   return 1; /* within bounds */
}



/* 
 * public objective function 
 * - parameters are on log scale
 * - if parameters are out of bounds function returns large value
 */
void objFuncLogScale(
 double *params,
 int *expo,
 double *result
)
{
  
 if( checkLogParams(params) == 0 ) /* out of bounds */
 {
   *result = LARGENUM;
 }
 else
 {
   nprsObjectiveFunction(params, expo, result);
 }
}


/*============================================================================*/


double ageAlongEdges(int node) { /* tree must be clock-like */
 int i;

 if(node>=0) return(0.);

 for(i=0;i<tree_nedges;i++)
  if(tree_lowerNodes[i]==node)
   return(tree_edgeLengths[i]+ageAlongEdges(tree_upperNodes[i])); 

 return(0.);
}

/* compute set of parameters for a given clock-like tree */
void getExternalParams(double *result) {
 int i;

 for(i=0;i<tree_nedges;i++) {
  if(tree_lowerNodes[i]<0&&tree_upperNodes[i]<0)
   result[index[-tree_upperNodes[i]]]
    =-log(ageAlongEdges(tree_lowerNodes[i])/ageAlongEdges(tree_upperNodes[i]));
 }

}

/*============================================================================*/
