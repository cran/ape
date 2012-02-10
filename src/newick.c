/* newick.c    2012-02-09 */

/* Copyright 2007-2008 Vincent Lefort */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

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
	      error("error reading subtree");
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
	      error("node %s has too many (>3) children.", centerNode->label);
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
      error("generated tree does not start with '('");
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
	  error("generated tree has too many right parentheses");
	}
    }
  centerNode = decodeNewickSubtree (str, T, &uCount);
  snprintf (centerNode->label, MAX_LABEL_LENGTH, "%s", rootLabel); /* added "%s" following Jos Kafer's suggestion (2010-11-23) */
  T->root = centerNode;
  return(T);
}
