/* BIONJ.c    2012-02-09 */

/* Copyright 2007-2008 Olivier Gascuel, Hoa Sien Cuong,
   R port by Vincent Lefort, bionj() below modified by Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                         BIONJ program                                     ;
;                                                                           ;
;                         Olivier Gascuel                                   ;
;                                                                           ;
;                         GERAD - Montreal- Canada                          ;
;                         olivierg@crt.umontreal.ca                         ;
;                                                                           ;
;                         LIRMM - Montpellier- France                       ;
;                         gascuel@lirmm.fr                                  ;
;                                                                           ;
;                         UNIX version, written in C                        ;
;                         by Hoa Sien Cuong (Univ. Montreal)                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

#include "me.h"

void Initialize(float **delta, double *X, char **labels, int n, POINTERS *trees);
void Print_outputChar(int i, POINTERS *trees, char *output);
void bionj(double *X, int *N, char **labels, int *edge1, int *edge2, double *el, char **tl);
int Symmetrize(float **delta, int n);
void Concatenate(char chain1[MAX_LABEL_LENGTH], int ind, POINTERS *trees, int post);
float Distance(int i, int j, float **delta);
float Variance(int i, int j, float **delta);
int Emptied(int i, float **delta);
float Sum_S(int i, float **delta);
void Compute_sums_Sx(float **delta, int n);
void Best_pair(float **delta, int r, int *a, int *b, int n);
float Finish_branch_length(int i, int j, int k, float **delta);
void FinishStr (float **delta, int n, POINTERS *trees, char *StrTree);
float Agglomerative_criterion(int i, int j, float **delta, int r);
float Branch_length(int a, int b, float **delta, int r);
float Reduction4(int a, float la, int b, float lb, int i, float lamda, float **delta);
float Reduction10(int a, int b, int i, float lamda, float vab, float **delta);
float Lamda(int a, int b, float vab, float **delta, int n, int r);

/* void printMat(float **delta, int n); */

/*;;;;;;;;;;;  INPUT, OUTPUT, INITIALIZATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                           ;
;                                                                           ;
;              The delta matrix is read from the input-file.                ;
;              It is recommended to put it and the executable in            ;
;              a special directory. The input-file and output-file          ;
;              can be given as arguments to the executable by               ;
;              typing them after the executable (Bionj input-file           ;
;              output-file) or by typing them when asked by the             ;
;              program. The input-file has to be formated according         ;
;              the PHYLIP standard. The output file is formated             ;
;              according to the NEWWICK standard.                           ;
;                                                                           ;
;              The lower-half of the delta matrix is occupied by            ;
;              dissimilarities. The upper-half of the matrix is             ;
;              occupied by variances. The first column                      ;
;              is initialized as 0; during the algorithm some               ;
;              indices are no more used, and the corresponding              ;
;              positions in the first column are set to 1.                  ;
;                                                                           ;
;              This delta matix is made symmetrical using the rule:         ;
;              Dij = Dji <- (Dij + Dji)/2. The diagonal is set to 0;        ;
;              during the further steps of the algorithm, it is used        ;
;              to store the sums Sx.                                        ;
;                                                                           ;
;              A second array, trees, is used to store taxon names.         ;
;              During the further steps of the algoritm, some               ;
;              positions in this array are emptied while the others         ;
;              are used to store subtrees.                                  ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


/*;;;;;;;;;;;;;;;;;;;;;;;;;; Initialize        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function reads an input file and return the            ;
;               delta matrix and trees: the list of taxa.                   ;
;                                                                           ;
; input       :                                                             ;
;              float **delta : delta matrix                                 ;
;              FILE *input    : pointer to input file                       ;
;              int n          : number of taxa                              ;
;              char **trees   : list of taxa                                ;
;                                                                           ;
; return value:                                                             ;
;              float **delta : delta matrix                                 ;
;              char *trees    : list of taxa                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

//void Initialize(float **delta, FILE *input, int n, POINTERS *trees)
void Initialize(float **delta, double *X, char **labels, int n, POINTERS *trees)
{
  int lig;                                          /* matrix line       */
  int col;                                          /* matrix column     */
//  float distance;
  //char name_taxon[LEN];                             /* taxon name        */
  char name_taxon[MAX_LABEL_LENGTH];
  char format[MAX_DIGITS];
  WORD *name;

  snprintf (format, MAX_DIGITS, "%%%ds", MAX_LABEL_LENGTH);

  for(lig=1; lig <= n; lig++)
    {
      //fscanf(input,"%s",name_taxon);                  /* read taxon name   */
      //fscanf (input, format, name_taxon);             /* read taxon name   */
      strncpy (name_taxon, labels[lig-1], MAX_LABEL_LENGTH);
      name=(WORD *)calloc(1,sizeof(WORD));            /* taxon name is     */
      if(name == NULL)                                /* put in trees      */
	{
	  error("out of memory");
	}
      else
	{
	  strncpy (name->name, name_taxon, MAX_LABEL_LENGTH);
	  name->suiv=NULL;
	  trees[lig].head=name;
	  trees[lig].tail=name;
	  for(col=lig; col <= n; col++)
	    {
	      //fscanf(input,"%f",&distance);             /* read the distance  */
//	      &distance = X[XINDEX(lig,col)];
	      delta[col][lig]=X[XINDEX(lig,col)];
	      delta[lig][col]=X[XINDEX(lig,col)];
	      if (lig==col)
	        delta[lig][col]=0;
	    }
	}
    }
  return;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Print_output;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                                                                           ;
; Description : This function prints out the tree in the output file.       ;
;                                                                           ;
; input       :                                                             ;
;              POINTERS *trees : pointer to the subtrees.                   ;
;              int i          : indicate the subtree i to be printed.       ;
:              FILE *output   : pointer to the output file.                 ;
;                                                                           ;
; return value: The phylogenetic tree in the output file.                   ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Print_outputChar(int i, POINTERS *trees, char *output)
{
  WORD *parcour;
  parcour=trees[i].head;
  while (parcour != NULL && (strlen (output) + strlen (parcour->name) < MAX_INPUT_SIZE))
    {
      output = strncat (output, parcour->name, strlen (parcour->name));
      parcour=parcour->suiv;
    }
  return;
}

//tree *bionj (FILE *input, boolean isNJ)
void bionj(double *X, int *N, char **labels,
	   int *edge1, int *edge2, double *el, char **tl)
{
  POINTERS *trees;            /* list of subtrees            */
  tree *T;                    /* the returned tree           */
  char *chain1;               /* stringized branch-length    */
  char *str;                  /* the string containing final tree */
  int *a, *b;                 /* pair to be agglomerated     */
  float **delta;              /* delta matrix                */
  float la;                   /* first taxon branch-length   */
  float lb;                   /* second taxon branch-length  */
  float vab;                  /* variance of Dab             */
  float lamda = 0.5;
  int i;
//  int ok;
  int r;                      /* number of subtrees          */
  int n;                      /* number of taxa              */
  int x, y;
//  double t;
  a=(int*)calloc(1,sizeof(int));
  b=(int*)calloc(1,sizeof(int));
  chain1=(char *)R_alloc(MAX_LABEL_LENGTH, sizeof(char));
  str = (char *)R_alloc(MAX_INPUT_SIZE, sizeof(char));
  /* added by EP */
  if (strlen(chain1))
    strncpy(chain1, "", strlen(chain1));
  if (strlen(str))
    strncpy(str, "", strlen(str));
  /* end */

//  fscanf(input,"%d",&n);
  n = *N;
  /*      Create the delta matrix     */
  delta=(float **)calloc(n+1,sizeof(float*));
  for(i=1; i<= n; i++)
    {
      delta[i]=(float *)calloc(n+1, sizeof(float));
      if(delta[i] == NULL)
	{
	  error("out of memory");
	}
    }
  trees=(POINTERS *)calloc(n+1,sizeof(POINTERS));
  if(trees == NULL)
    {
      error("out of memory");
    }
  /*   initialise and symmetrize the running delta matrix    */
  r=n;
  *a=0;
  *b=0;
  Initialize(delta, X, labels, n, trees);
//  ok=Symmetrize(delta, n);

//  if(!ok)
//   Rprintf("\n The matrix is not symmetric.\n ");
  while (r > 3)                             /* until r=3                 */
      {
      	Compute_sums_Sx(delta, n);             /* compute the sum Sx       */
	Best_pair(delta, r, a, b, n);          /* find the best pair by    */
	vab=Variance(*a, *b, delta);           /* minimizing (1)           */
	la=Branch_length(*a, *b, delta, r);    /* compute branch-lengths   */
	lb=Branch_length(*b, *a, delta, r);    /* using formula (2)        */
//	if (!isNJ)
	  lamda=Lamda(*a, *b, vab, delta, n, r); /* compute lambda* using (9)*/
	for(i=1; i <= n; i++)
	  {
	    if(!Emptied(i,delta) && (i != *a) && (i != *b))
	      {
		if(*a > i)
		  {
		    x=*a;
		    y=i;
		  }
		else
		  {
		    x=i;
		    y=*a;                           /* apply reduction formulae */
		  }                                 /* 4 and 10 to delta        */
		delta[x][y]=Reduction4(*a, la, *b, lb, i, lamda, delta);
		delta[y][x]=Reduction10(*a, *b, i, lamda, vab, delta);
	      }
	  }
	strncpy(chain1,"",1);                  /* agglomerate the subtrees */
	strncat(chain1,"(",1);                 /* a and b together with the*/
	Concatenate(chain1, *a, trees, 0);     /* branch-lengths according */
	strncpy(chain1,"",1);                  /* to the NEWWICK format    */
	strncat(chain1,":",1);
	snprintf(chain1+strlen(chain1),MAX_LABEL_LENGTH,"%f",la);

	strncat(chain1,",",1);
	Concatenate(chain1,*a, trees, 1);
	trees[*a].tail->suiv=trees[*b].head;
	trees[*a].tail=trees[*b].tail;
	strncpy(chain1,"",1);
	strncat(chain1,":",1);
	snprintf(chain1+strlen(chain1),MAX_LABEL_LENGTH,"%f",lb);

	strncat(chain1,")",1);
	Concatenate(chain1, *a, trees, 1);
	delta[*b][0]=1.0;                     /* make the b line empty     */
	trees[*b].head=NULL;
	trees[*b].tail=NULL;
	r=r-1;
      }

  FinishStr (delta, n, trees, str);   /* compute the branch-lengths*/
                                      /* of the last three subtrees*/
				      /* and print the tree in the */
				      /* output-file               */
  T = readNewickString (str, n);
  T = detrifurcate(T);
//  compareSets(T,species);
  partitionSizes(T);

  tree2phylo(T, edge1, edge2, el, tl, n); /* by EP */

  for(i=1; i<=n; i++)
  {
      delta[i][0]=0.0;
      trees[i].head=NULL;
      trees[i].tail=NULL;
  }
  free(delta);
  free(trees);
  /* free (str); */
  /* free (chain1); */
  free(a);
  free(b);
  freeTree(T);
  T = NULL;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                             Utilities                                     ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Symmetrize  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function verifies if the delta matrix is symmetric;    ;
;               if not the matrix is made symmetric.                        ;
;                                                                           ;
; input       :                                                             ;
;              float **delta : delta matrix                                 ;
;              int n          : number of taxa                              ;
;                                                                           ;
; return value:                                                             ;
;              int symmetric  : indicate if the matrix has been made        ;
;                               symmetric or not                            ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

int Symmetrize(float **delta, int n)
{
  int lig;                                         /* matrix line        */
  int col;                                         /* matrix column      */
  float value;                                     /* symmetrized value  */
  int symmetric;

  symmetric=1;
  for(lig=1; lig  <=  n; lig++)
    {
      for(col=1; col < lig; col++)
	{
	  if(delta[lig][col] != delta[col][lig])
	    {
	      value= (delta[lig][col]+delta[col][lig])/2;
	      delta[lig][col]=value;
	      delta[col][lig]=value;
	      symmetric=0;
	    }
        }
    }
  return(symmetric);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;; Concatenate ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                                                                           ;
; Description : This function concatenates a string to another.             ;
;                                                                           ;
; input       :                                                             ;
;      char *chain1    : the string to be concatenated.                     ;
;      int ind         : indicate the subtree to which concatenate the      ;
;                        string                                             ;
;      POINTERS *trees  : pointer to subtrees.                              ;
;      int post        : position to which concatenate (front (0) or        ;
;                        end (1))                                           ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

//void Concatenate(char chain1[LEN], int ind, POINTERS *trees, int post)
void Concatenate(char chain1[MAX_LABEL_LENGTH], int ind, POINTERS *trees, int post)
{
  WORD *bran;

  bran=(WORD *)calloc(1,sizeof(WORD));
  if(bran == NULL)
    {
      error("out of memory");
    }
  else
    {
      strncpy(bran->name,chain1,MAX_LABEL_LENGTH);
      bran->suiv=NULL;
    }
  if(post == 0)
    {
      bran->suiv=trees[ind].head;
      trees[ind].head=bran;
    }
  else
    {
      trees[ind].tail->suiv=bran;
      trees[ind].tail=trees[ind].tail->suiv;
    }
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Distance;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function retrieve ant return de distance between taxa  ;
;               i and j from the delta matrix.                              ;
;                                                                           ;
; input       :                                                             ;
;               int i          : taxon i                                    ;
;               int j          : taxon j                                    ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
; return value:                                                             ;
;               float distance : dissimilarity between the two taxa         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Distance(int i, int j, float **delta)
{
  if(i > j)
    return(delta[i][j]);
  else
    return(delta[j][i]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Variance;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function retrieve and return the variance of the       ;
;               distance between i and j, from the delta matrix.            ;
;                                                                           ;
; input       :                                                             ;
;               int i           : taxon i                                   ;
;               int j           : taxon j                                   ;
;               float **delta  : the delta matrix                           ;
;                                                                           ;
; return value:                                                             ;
;               float distance : the variance of  Dij                       ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Variance(int i, int j, float **delta)
{
  if(i > j)
    return(delta[j][i]);
  else
    return(delta[i][j]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Emptied ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function verifie if a line is emptied or not.          ;
;                                                                           ;
; input       :                                                             ;
;               int i          : subtree (or line) i                        ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
; return value:                                                             ;
;               0              : if not emptied.                            ;
;               1              : if emptied.                                ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

int Emptied(int i, float **delta)      /* test if the ith line is emptied */
{
  return((int)delta[i][0]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Sum_S;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function retrieves the sum Sx from the diagonal       ;
;                of the delta matrix.                                       ;
;                                                                           ;
;  input       :                                                            ;
;               int i          : subtree i                                  ;
;               float **delta : the delta matrix                            ;
;                                                                           ;
;  return value:                                                            ;
;                float delta[i][i] : sum Si                                 ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Sum_S(int i, float **delta)          /* get sum Si form the diagonal */
{
  return(delta[i][i]);
}


/*;;;;;;;;;;;;;;;;;;;;;;;Compute_sums_Sx;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
; Description : This function computes the sums Sx and store them in the    ;
;               diagonal the delta matrix.                                  ;
;                                                                           ;
; input       :                                                             ;
;     	         float **delta : the delta matrix.                      ;
;     	         int n          : the number of taxa                    ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Compute_sums_Sx(float **delta, int n)
{
  float sum;
  int i;
  int j;

  for(i= 1; i <= n ; i++)
    {
      if(!Emptied(i,delta))
	{
	  sum=0;
	  for(j=1; j <=n; j++)
	    {
	      if(i != j && !Emptied(j,delta))           /* compute the sum Si */
		sum=sum + Distance(i,j,delta);
	    }
	}
      delta[i][i]=sum;                           /* store the sum Si in */
    }                                               /* delta’s diagonal    */
}


/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Best_pair;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description : This function finds the best pair to be agglomerated by    ;
;                minimizing the agglomerative criterion (1).                ;
;                                                                           ;
;  input       :                                                            ;
;                float **delta : the delta matrix                           ;
;                int r          : number of subtrees                        ;
;                int *a         : contain the first taxon of the pair       ;
;                int *b         : contain the second taxon of the pair      ;
;                int n          : number of taxa                            ;
;                                                                           ;
;  return value:                                                            ;
;                int *a         : the first taxon of the pair               ;
;                int *b         : the second taxon of the pair              ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

void Best_pair(float **delta, int r, int *a, int *b, int n)
{
  float Qxy;                         /* value of the criterion calculated*/
  int x,y;                           /* the pair which is tested         */
  float Qmin;                        /* current minimun of the criterion */

  Qmin=1.0e300;
  for(x=1; x <= n; x++)
    {
      if(!Emptied(x,delta))
        {
	  for(y=1; y < x; y++)
	    {
	      if(!Emptied(y,delta))
		{
		  Qxy=Agglomerative_criterion(x,y,delta,r);
		  if(Qxy < Qmin-0.000001)
		    {
		      Qmin=Qxy;
		      *a=x;
		      *b=y;
		    }
		}
	    }
        }
    }
}


/*;;;;;;;;;;;;;;;;;;;;;;Finish_branch_length;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;  Description :  Compute the length of the branch attached                 ;
;                 to the subtree i, during the final step                   ;
;                                                                           ;
;  input       :                                                            ;
;                int i          : position of subtree i                     ;
;                int j          : position of subtree j                     ;
;                int k          : position of subtree k                     ;
;                float **delta :                                            ;
;                                                                           ;
;  return value:                                                            ;
;                float length  : The length of the branch                   ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/

float Finish_branch_length(int i, int j, int k, float **delta)
{
  float length;
  length=0.5*(Distance(i,j,delta) + Distance(i,k,delta)
	      -Distance(j,k,delta));
  return(length);
}

void FinishStr (float **delta, int n, POINTERS *trees, char *StrTree)
{
  int l=1;
  int i=0;
  float length;
  char *tmp;
  WORD *bidon;
  WORD *ele;
  int last[3];                            /* the last three subtrees     */

  while(l <= n)
    {                                     /* find the last tree subtree  */
      if(!Emptied(l, delta))
	{
	  last[i]=l;
	  i++;
	}
      l++;
    }
  tmp = (char*) calloc (12, sizeof(char));
  StrTree[0]='(';

  length=Finish_branch_length(last[0],last[1],last[2],delta);
  Print_outputChar (last[0], trees, StrTree);
  snprintf (tmp, 12, "%f,", length);
  if (strlen (StrTree) + strlen (tmp) < MAX_INPUT_SIZE-1) {
    strncat (StrTree, ":", 1);
    strncat (StrTree, tmp, strlen (tmp));
  }

  length=Finish_branch_length(last[1],last[0],last[2],delta);
  Print_outputChar (last[1], trees, StrTree);
  snprintf (tmp, 12, "%f,", length);
  if (strlen (StrTree) + strlen (tmp) < MAX_INPUT_SIZE-1) {
    strncat (StrTree, ":", 1);
    strncat (StrTree, tmp, strlen (tmp));
  }

  length=Finish_branch_length(last[2],last[1],last[0],delta);
  Print_outputChar (last[2], trees, StrTree);
  snprintf (tmp, 12, "%f", length);
  if (strlen (StrTree) + strlen (tmp) < MAX_INPUT_SIZE-1) {
    strncat (StrTree, ":", 1);
    strncat (StrTree, tmp, strlen (tmp));
  }

  if (strlen (StrTree) < MAX_INPUT_SIZE-3)
    strncat (StrTree, ");", 3);

  free (tmp);
  for(i=0; i < 3; i++)
    {
      bidon=trees[last[i]].head;
      ele=bidon;
      while(bidon!=NULL)
	{
	  ele=ele->suiv;
	  free(bidon);
	  bidon=ele;
	}
    }
  return;
}

/*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*\
;                                                                           ;
;                          Formulae                                         ;
;                                                                           ;
\*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;*/


float Agglomerative_criterion(int i, int j, float **delta, int r)
{
  float Qij;
  Qij=(r-2)*Distance(i,j,delta)                           /* Formula (1) */
    -Sum_S(i,delta)
    -Sum_S(j,delta);

  return(Qij);
}


float Branch_length(int a, int b, float **delta, int r)
{
  float length;
  length=0.5*(Distance(a,b,delta)                         /* Formula (2) */
	      +(Sum_S(a,delta)
		-Sum_S(b,delta))/(r-2));
  return(length);
}


float Reduction4(int a, float la, int b, float lb, int i, float lamda, float **delta)
{
  float Dui;
  Dui=lamda*(Distance(a,i,delta)-la)
    +(1-lamda)*(Distance(b,i,delta)-lb);                /* Formula (4) */
  return(Dui);
}


float Reduction10(int a, int b, int i, float lamda, float vab, float **delta)
{
  float Vci;
  Vci=lamda*Variance(a,i,delta)+(1-lamda)*Variance(b,i,delta)
    -lamda*(1-lamda)*vab;                              /*Formula (10)  */
  return(Vci);
}

float Lamda(int a, int b, float vab, float **delta, int n, int r)
{
  float lamda=0.0;
  int i;

  if(vab==0.0)
    lamda=0.5;
  else
    {
      for(i=1; i <= n ; i++)
	{
          if(a != i && b != i && !Emptied(i,delta))
            lamda=lamda + (Variance(b,i,delta) - Variance(a,i,delta));
	}
      lamda=0.5 + lamda/(2*(r-2)*vab);
    }                                       /* Formula (9) and the  */
  if(lamda > 1.0)                           /* constraint that lamda*/
    lamda = 1.0;                            /* belongs to [0,1]     */
  if(lamda < 0.0)
    lamda=0.0;
  return(lamda);
}

/*
void printMat(float **delta, int n)
{
  int i, j;
  for (i=1; i<=n; i++) {
    for (j=1; j<=n; j++)
      Rprintf ("%f ", delta[i][j]);
    Rprintf("\n");
  }
  Rprintf("\n");
  return;
}
*/
