
/* 
 *  treefunc.c
 *
 * (c) 2003  Gangolf Jobb (http://www.treefinder.de)
 *
 *  Various data structures and methods for manipulating
 *  and traversing trees
 *  (e.g., to convert between "phylo" and "hclust" objects
 *  and to classify genes)
 *  
 *
 *  This code may be distributed under the GNU GPL
 */


#include <stdlib.h>
#include <string.h>
#include <math.h>

/*============================================================================*/

#define LABLEN 128

typedef struct tree {
 char label[LABLEN];
 struct tree *branches,*next;
 double length;
 int mark;
} Tree;

/*....*/

Tree *NewTree(void) {
 Tree *b;

 b=malloc(sizeof(Tree)); if(!b) return(NULL);

 *b->label='\0'; b->branches=NULL; b->next=NULL; b->length=0.0;

 return(b);
}

/*....*/

void FreeTree(Tree *t) {
 Tree *b;

 while(t->branches) {b=t->branches; t->branches=b->next; FreeTree(b);}

 free(t);

}

/*============================================================================*/

static Tree *current=NULL;

static int tip_index;
static int edge_index; 
static int node_index;

enum {OK=0,ERROR}; /* one may invent more descriptive error codes */
static int error=OK;

/*============================================================================*/

Tree *buildTreeFromPhylo_(
 int node,
 int *lowerNodes,
 int *upperNodes,
 double *edgeLengths,
 int nedges,
 char **tipLabels,
 int ntips
) {
 Tree *t,*b,**bb;
 int i,j,n;

 t=NewTree();

 bb=&t->branches; n=0;
 for(i=0;i<nedges;i++) { if(lowerNodes[i]!=node) continue;
  j=upperNodes[i];                                             if(j==0) {error=ERROR; goto err;}
  if(j>0) {                                                    if(j>ntips) {error=ERROR; goto err;}
   b=NewTree(); strcpy(b->label,tipLabels[j-1]);
  } else {                                                     if(-j>nedges) {error=ERROR; goto err;}
   b=
    buildTreeFromPhylo_(j,lowerNodes,upperNodes,edgeLengths,nedges,tipLabels,ntips);
  }
  b->length=edgeLengths[i];
  *bb=b; bb=&b->next; n++;
 }                                                             if(n<2) {error=ERROR; goto err;}
                                                               err:
 *bb=NULL;

 return(t);
}

/*....*/

void buildTreeFromPhylo(
 int *lowerNodes,
 int *upperNodes,
 double *edgeLengths,
 int *nedges,
 char **tipLabels,
 int *ntips,
 int *result
) {

 error=OK;

 if(current) {FreeTree(current); current=NULL;}

 if(*nedges<2||*ntips<2) {error=ERROR; *result=error; return;}

 current=buildTreeFromPhylo_(-1,lowerNodes,upperNodes,edgeLengths,*nedges,tipLabels,*ntips);

 if(error&&current) {FreeTree(current); current=NULL;}

 *result=error;

}

/*============================================================================*/

void destroyTree(int *result) {

 error=OK;

 if(current) {FreeTree(current); current=NULL;}

 *result=error;
 
}

/*============================================================================*/

void getError(int *result) {

 *result=error;

 error=OK;

}

/*============================================================================*/

int nTips_(Tree *t) {
 Tree *b;
 int n;

 if(!t->branches) return(1);

 n=0; for(b=t->branches;b;b=b->next) n+=nTips_(b);

 return(n);
}

/*....*/

void nTips(int *result) {

 error=OK;

 if(!current) {error=ERROR; *result=0; return;}

 *result=nTips_(current);
 
}

/*============================================================================*/

int nNodes_(Tree *t) {
 Tree *b;
 int n;

 if(!t->branches) return(1);

 n=1; for(b=t->branches;b;b=b->next) n+=nNodes_(b);

 return(n);
}

/*....*/

void nNodes(int *result) {

 error=OK;

 if(!current) {error=ERROR; *result=0; return;}

 *result=nNodes_(current);

}

/*....*/

void nEdges(int *result) {

 error=OK;

 if(!current) {error=ERROR; *result=0; return;}

 *result=nNodes_(current)-1;
 
}

/*============================================================================*/

double markClasses_(Tree *t) {
 Tree *b; double destinct,sum;

 /* all tips above a marked ( == 1) node belong to the same class */

 if(!t->branches) {t->mark=1; return(t->length);}

 sum=0.; for(b=t->branches;b;b=b->next) sum+=markClasses_(b);

 destinct=nTips_(t)*(t->length); /* (t->length) == 0. at root */

 if(destinct>sum) { t->mark=1;  return(destinct); }

                    t->mark=0;  return(sum);
 
}

/*....*/

void getMisawaTajima__(Tree *t,int ignore,int *result) { /* maps tips to marked classes */
 Tree *b;

 if(t->mark&&!ignore) {node_index++; ignore=1;} /* marked nodes above a marked node will be ignored */

 if(!t->branches) {result[tip_index++]=node_index; return;}

 for(b=t->branches;b;b=b->next) getMisawaTajima__(b,ignore,result);

}

/*....*/

void getMisawaTajima_(Tree *t,int *result) {

 markClasses_(t);

 tip_index=0; node_index=0;

 getMisawaTajima__(t,0,result);

}

/*....*/

void getMisawaTajima(int *result) {

 error=OK;

 if(!current) {error=ERROR; return;}

 getMisawaTajima_(current,result);
 
}
