/*#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "graph.h"
#include "newick.h"
#include "main.h"
*/
#include "me.h"

int nodeCount;
int edgeCount;

int whiteSpace(char c)
{
  if ((' ' == c) || ('\t' == c) || ('\n' == c))
    return(1);
  else return(0);
}

int leaf(node *v)
{
  int count = 0;
  if (NULL != v->parentEdge)
    count++;
  if (NULL != v->leftEdge)
    count++;
  if (NULL != v->rightEdge)
    count++;
  if (NULL != v->middleEdge)
    count++;
  if (count > 1)
    return(0);
  return(1);
}

/*decodeNewickSubtree is used to turn a string of the form
  "(v1:d1,v2:d2,(subtree) v3:d3....vk:dk) subroot:d," into a subtree
  rooted at subrooted, with corresponding subtrees and leaves at v1
  through vk.  It is called recursively on subtrees*/

node *decodeNewickSubtree(char *treeString, tree *T, int *uCount)
{
  node *thisNode;
  node *centerNode;
  double thisWeight;
  edge *thisEdge;
//  char label[MAX_LABEL_LENGTH];
  char stringWeight[MAX_LABEL_LENGTH];
  int state;
  int i = 0;
  int j;
  int left,right;
  int parcount;
//  snprintf(label,14,"Default_Label\0");
  left = right = 0;
  parcount = 0;
  state = ReadOpenParenthesis;
  if('(' == treeString[0])
    parcount++;
  //centerNode = makeNode(label,NULL,nodeCount++);
  centerNode = makeNode("",NULL,nodeCount++);
  T->size++;
  while(parcount > 0)
    {
      while(whiteSpace(treeString[i]))
	i++;
      switch(state)
	{
	case(ReadOpenParenthesis):
	  if('(' != treeString[0])
	    {
	      Rprintf("Error reading subtree.\n");
	      exit(0);
	    }
	  i++;
	  state = ReadSubTree;
	  break;
	case(ReadSubTree):
	  if('(' == treeString[i])  /*if treeString[i] is a left parenthesis,
				      we scan down the string until we find its partner.
				      the relative number of '('s vs. ')'s is counted
				      by the variable parcount*/
	    {
	      left = i++;
	      parcount++;
	      while(parcount > 1)
		{
		  while (('(' != treeString[i]) && (')' != treeString[i]))
		    i++;  /*skip over characters which are not parentheses*/
		  if('(' == treeString[i])
		    parcount++;
		  else if (')' == treeString[i])
		    parcount--;
		  i++;
		}  /*end while */
	      right = i;  /*at this point, the subtree string goes from
			    treeString[left] to treeString[right - 1]*/
	      thisNode = decodeNewickSubtree(treeString + left,T,uCount);  /*note that this
								      step will put
								      thisNode in T*/
	      i = right;  /*having created the node for the subtree, we move
			    to assigning the label for the new node.
			    treeString[right] will be the start of this label */
	    } /* end if ('(' == treeString[i]) condition */
	  else
	    {
	      //thisNode = makeNode(label,NULL,nodeCount++);
	      thisNode = makeNode("",NULL,nodeCount++);
	      T->size++;
	    }
	  state = ReadLabel;
	  break;
	case(ReadLabel):
	  left = i;  /*recall "left" is the left marker for the substring, "right" the right*/
	  if (':' == treeString[i])   /*this means an internal node?*/
	    {
	      //sprintf(thisNode->label,"I%d",(*uCount)++);
	      //snprintf(thisNode->label,MAX_LABEL_LENGTH,"I%d",(*uCount)++);
	      (*uCount)++;
	      right = i;
	    }
	  else
	    {
	      while((':' != treeString[i]) && (',' != treeString[i]) && (')' != treeString[i]))
		i++;
	      right = i;
	      j = 0;
	      for(i = left; i < right;i++)
		if(!(whiteSpace(treeString[i])))
		  thisNode->label[j++] = treeString[i];
	      thisNode->label[j] = '\0';
	    }
	  if(':' == treeString[right])
	    state = ReadWeight;
	  else
	    {
	      state = AddEdge;
	      thisWeight = 0.0;
	    }
	  i = right + 1;
	  break;
	case(ReadWeight):
	  left = i;
	  while
	    ((',' != treeString[i]) && (')' != treeString[i]))
	    i++;
	  right = i;
	  j = 0;
	  for(i = left; i < right; i++)
	    stringWeight[j++] = treeString[i];
	  stringWeight[j] = '\0';
	  thisWeight = atof(stringWeight);
	  state=AddEdge;
	  break;
	case(AddEdge):
	  //thisEdge = makeEdge(label,centerNode,thisNode,thisWeight);
	  thisEdge = makeEdge("",centerNode,thisNode,thisWeight);
	  thisNode->parentEdge = thisEdge;
	  if (NULL == centerNode->leftEdge)
	    centerNode->leftEdge = thisEdge;
	  else if (NULL == centerNode->rightEdge)
	    centerNode->rightEdge = thisEdge;
	  else if (NULL == centerNode->middleEdge)
	    centerNode->middleEdge = thisEdge;
	  else
	    {
	      Rprintf("Error: node %s has too many (>3) children.\n",centerNode->label);
	      exit(0);
	    }
	  //sprintf(thisEdge->label,"E%d",edgeCount++);
	  //snprintf(thisEdge->label,MAX_LABEL_LENGTH,"E%d",edgeCount++);
	  edgeCount++;
	  i = right + 1;
	  if (',' == treeString[right])
	    state = ReadSubTree;
	  else
	    parcount--;
	  break;
	}
    }
  return(centerNode);
}

tree *readNewickString (char *str, int numLeaves)
{
  tree *T;
  node *centerNode;
  int i = 0;
  int j = 0;
  int inputLength;
  int uCount = 0;
  int parCount = 0;
  char rootLabel[MAX_LABEL_LENGTH];
  nodeCount = edgeCount = 0;

  T = newTree();

  if ('(' != str[0])
    {
      Rprintf("Error reading generated tree - does not start with '('.\n");
      exit(0);
    }
  inputLength = strlen (str)+1;
  for(i = 0; i < inputLength; i++)
    {
      if ('(' == str[i])
	parCount++;
      else if (')' == str[i])
	parCount--;
      if (parCount > 0)
	;
      else if (0 == parCount)
	{
	  i++;
/*	  if(';' == str[i])
	    sprintf(rootLabel,"URoot");
	  else
	    {*/
	      while(';' != str[i])
	        if (!(whiteSpace (str[i++])))
	          rootLabel[j++] = str[i-1];  /*be careful here */
	        rootLabel[j] = '\0';
//	    }
	  i = inputLength;
	}
      else if (parCount < 0)
	{
	  Rprintf("Error reading generated tree. Too many right parentheses.\n");
	  exit(0);
	}
    }
  centerNode = decodeNewickSubtree (str, T, &uCount);
  snprintf (centerNode->label, MAX_LABEL_LENGTH, rootLabel);
  T->root = centerNode;
  return (T);
}

tree *loadNewickTree(FILE *ifile, int numLeaves)
{
//  char label[] = "EmptyEdge";
  tree *T;
  node *centerNode;
  int i = 0;
  int j = 0;
  int inputLength;
  int uCount = 0;
  int parCount = 0;
  char c;
  int Comment;
  char *nextString;
  char rootLabel[MAX_LABEL_LENGTH];
  nodeCount = edgeCount = 0;
  T = newTree();
  nextString = (char *) malloc(numLeaves*INPUT_SIZE*sizeof(char));
  if (NULL == nextString)
    nextString = (char *) malloc(MAX_INPUT_SIZE*sizeof(char));
  Comment = 0;
  while(1 == fscanf(ifile,"%c",&c))
    {
      if('[' == c)
	Comment = 1;
      else if (']' == c)
	Comment = 0;
      else if (!(Comment))
	{
	  if(whiteSpace(c))
	    {
	      if (i > 0)
		nextString[i++] = ' ';
	    }
	  else  /*note that this else goes with if(whiteSpace(c))*/
	    nextString[i++] = c;
	  if (';' == c)
	    break;
	}
    }
  if ('(' != nextString[0])
    {
      fprintf(stderr,"Error reading input file - does not start with '('.\n");
      exit(EXIT_FAILURE);
    }
  inputLength = i;
  for(i = 0; i < inputLength;i++)
    {
      if ('(' == nextString[i])
	parCount++;
      else if (')' == nextString[i])
	parCount--;
      if (parCount > 0)
	;
      else if (0 == parCount)
	{
	  i++;
/*	  if(';' == nextString[i])
	    sprintf(rootLabel,"URoot");
	  else
	    {*/
	      while(';' != nextString[i])
		if(!(whiteSpace(nextString[i++])))
		  rootLabel[j++] = nextString[i-1];  /*be careful here */
	      rootLabel[j] = '\0';
//	    }
	  i = inputLength;
	}
      else if (parCount < 0)
	{
	  fprintf(stderr,"Error reading tree input file.  Too many right parentheses.\n");
	  exit(EXIT_FAILURE);
	}
    }
  centerNode = decodeNewickSubtree(nextString,T,&uCount);
  snprintf(centerNode->label, MAX_LABEL_LENGTH, rootLabel);
  T->root = centerNode;
  free(nextString);
  return(T);
}

double GetSubTreeLength (tree *T, edge *e)
{
  double ret = 0;

  if ( (NULL != e) && (! leaf(e->head) )) {
    ret += GetSubTreeLength (T, e->head->leftEdge);
    ret += GetSubTreeLength (T, e->head->rightEdge);
  }
  ret += e->distance;
  return ret;
}

void NewickPrintSubtree(tree *T, edge *e, char *str)
{
  char *tmp;
  if (NULL == e)
    {
      Rprintf("Error with Newick Printing routine.\n");
      exit(0);
    }
  if(!(leaf(e->head)))
    {
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, "(", 1);
      NewickPrintSubtree(T,e->head->leftEdge,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ",", 1);
      NewickPrintSubtree(T,e->head->rightEdge,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ")", 1);
    }
  if (strlen (str) < MAX_INPUT_SIZE - strlen (e->head->label) -1)
    strncat (str, e->head->label, strlen (e->head->label));

  if (strlen (str) < MAX_INPUT_SIZE - 2)
    strncat (str, ":", 1);

  tmp = (char *)R_alloc(INPUT_SIZE, sizeof(char));
  /* added by EP */
  if (strlen(tmp))
    strncpy(tmp, "", strlen(tmp));
  /* end */
  snprintf (tmp, INPUT_SIZE, "%lf", e->distance);
  if (strlen (str) < MAX_INPUT_SIZE - strlen (tmp) -1)
    strncat (str, tmp, strlen (tmp));

  /* free (tmp); */
  return;
}

double GetBinaryTreeLength (tree *T)
{
  double ret = 0;
  edge *e, *f;
  node *rootchild;
  e = T->root->leftEdge;
  rootchild = e->head;

  f = rootchild->leftEdge;
  if (NULL != f)
    ret += GetSubTreeLength (T, f);
  f = rootchild->rightEdge;
  if (NULL != f)
    ret += GetSubTreeLength (T, f);
  ret += e->distance;
  return ret;
}

void NewickPrintBinaryTree(tree *T, char *str)
{
  edge *e, *f;
  node *rootchild;
  char *tmp;
  e = T->root->leftEdge;
  rootchild = e->head;
  if (strlen (str) < MAX_INPUT_SIZE -2)
    strncat (str, "(", 1);
  f = rootchild->leftEdge;
  if (NULL != f)
    {
      NewickPrintSubtree(T,f,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ",", 1);
    }
  f = rootchild->rightEdge;
  if (NULL != f)
    {
      NewickPrintSubtree(T,f,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ",", 1);
    }
  if (strlen (str) < MAX_INPUT_SIZE - strlen (T->root->label) -1)
    strncat (str, T->root->label, strlen (T->root->label));

  if (strlen (str) < MAX_INPUT_SIZE - 2)
    strncat (str, ":", 1);

  tmp = (char *)R_alloc(INPUT_SIZE, sizeof(char));
  /* added by EP */
  if (strlen(tmp))
    strncpy(tmp, "", strlen(tmp));
  /* end */
  snprintf (tmp, INPUT_SIZE, "%lf", e->distance);
  if (strlen (str) < MAX_INPUT_SIZE - strlen (tmp) -1)
    strncat (str, tmp, strlen (tmp));

  if (strlen (str) < MAX_INPUT_SIZE - 2)
    strncat (str, ")", 1);

  if (NULL != rootchild->label)
    if (strlen (str) < MAX_INPUT_SIZE - strlen (rootchild->label) -1)
      strncat (str, T->root->label, strlen (rootchild->label));

  if (strlen (str) < MAX_INPUT_SIZE - 3)
    strncat (str, ";\n", 2);

  /* free (tmp); */
  return;
}

double GetTrinaryTreeLength (tree *T)
{
  double ret = 0;
  edge *f;
  f = T->root->leftEdge;
  if (NULL != f)
    ret += GetSubTreeLength (T, f);
  f = T->root->rightEdge;
  if (NULL != f)
    ret += GetSubTreeLength (T, f);
  f = T->root->middleEdge;
  if (NULL != f)
    ret += GetSubTreeLength (T, f);

  return ret;
}

void NewickPrintTrinaryTree(tree *T, char *str)
{
  edge *f;
  f = T->root->leftEdge;
  if (strlen (str) < MAX_INPUT_SIZE -2)
    strncat (str, "(", 1);
  if (NULL != f)
    {
      NewickPrintSubtree(T,f,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ",", 1);
    }
  f = T->root->rightEdge;
  if (NULL != f)
    {
      NewickPrintSubtree(T,f,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ",", 1);
    }
  f = T->root->middleEdge;
  if (NULL != f)
    {
      NewickPrintSubtree(T,f,str);
      if (strlen (str) < MAX_INPUT_SIZE -2)
        strncat (str, ")", 1);
    }
  if (NULL != T->root->label)
    if (strlen (str) < MAX_INPUT_SIZE - strlen (T->root->label) -1)
      strncat (str, T->root->label, strlen (T->root->label));
  if (strlen (str) < MAX_INPUT_SIZE - 3)
    strncat (str, ";\n", 2);
  return;
}

void NewickPrintTreeStr(tree *T, char *str)
{
  if (leaf(T->root))
    NewickPrintBinaryTree(T,str);
  else
    NewickPrintTrinaryTree(T,str);
}

double GetTreeLength (tree *T)
{
  double ret = 0;
  if (leaf(T->root))
    ret = GetBinaryTreeLength (T);
  else
    ret = GetTrinaryTreeLength (T);
  return ret;
}
/*
void NewickPrintTree(tree *T, FILE *ofile)
{
  if (leaf(T->root))
    NewickPrintBinaryTree(T,ofile);
  else
    NewickPrintTrinaryTree(T,ofile);
}
*/
//edge *depthFirstTraverse(tree *T, edge *e);

