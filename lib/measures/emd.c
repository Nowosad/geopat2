/*
    emd.c

    Changes 12/15/2016: 
	adaptation to GeoPAT 2.0 (GRASS free)
    Author: Pawel Netzel, Space Informatics Lab

    Changes 11/05/2014: 
	cost matrix transfer from similarity library,
	error handling from GRASS,
	using dynamically allocated arrays.
    Author: Pawel Netzel, Space Informatics Lab

    Last update: 3/14/98

    An implementation of the Earth Movers Distance.
    Based of the solution for the Transportation problem as described in
    "Introduction to Mathematical Programming" by F. S. Hillier and 
    G. J. Lieberman, McGraw-Hill, 1990.

    Copyright (C) 1998 Yossi Rubner
    Computer Science Department, Stanford University
    E-Mail: rubner@cs.stanford.edu   URL: http://vision.stanford.edu/~rubner
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
#include <grass/glocale.h>
#include <grass/gis.h>
*/
#include "emd.h"
//#include "macro.h"
#include "geopat_compatibility.h"

#define DEBUG_LEVEL 0
/*
 DEBUG_LEVEL:
   0 = NO MESSAGES
   1 = PRINT THE NUMBER OF ITERATIONS AND THE FINAL RESULT
   2 = PRINT THE RESULT AFTER EVERY ITERATION
   3 = PRINT ALSO THE FLOW AFTER EVERY ITERATION
   4 = PRINT A LOT OF INFORMATION (PROBABLY USEFUL ONLY FOR THE AUTHOR)
*/


#define MAX_SIG_SIZE1  (MAX_SIG_SIZE+1)  /* FOR THE POSIBLE DUMMY FEATURE */

/* NEW TYPES DEFINITION */

/* node1_t IS USED FOR SINGLE-LINKED LISTS */
typedef struct node1_t {
  int i;
  double val;
  struct node1_t *Next;
} node1_t;

/* node1_t IS USED FOR DOUBLE-LINKED LISTS */
typedef struct node2_t {
  int i, j;
  double val;
  struct node2_t *NextC;               /* NEXT COLUMN */
  struct node2_t *NextR;               /* NEXT ROW */
} node2_t;


typedef struct {
/* GLOBAL VARIABLE DECLARATION */
	int _n1, _n2;                          /* SIGNATURES SIZES */
	double *_C;                            /* THE COST MATRIX */
	node2_t _X[MAX_SIG_SIZE1*2];           /* THE BASIC VARIABLES VECTOR */
/* VARIABLES TO HANDLE _X EFFICIENTLY */
	node2_t *_EndX, *_EnterX;
	char _IsX[MAX_SIG_SIZE1][MAX_SIG_SIZE1];
	node2_t *_RowsX[MAX_SIG_SIZE1], *_ColsX[MAX_SIG_SIZE1];
	double _maxW;
	double _maxC;
} process_variables;

/* DECLARATION OF FUNCTIONS */
float emd(int size, double *hist1, double *hist2,
	  double *Dist, double max_dist,
	  flow_t *Flow, int *FlowSize);
double init(int size, double hist1[], double hist2[], double Dist[], double max_dist, process_variables *glob);
void findBasicVariables(node1_t *U, node1_t *V, process_variables *glob);
int isOptimal(node1_t *U, node1_t *V, process_variables *glob);
int findLoop(node2_t **Loop, process_variables *glob);
void newSol(process_variables *glob);
void russel(double *S, double *D, process_variables *glob);
void addBasicVariable(int minI, int minJ, double *S, double *D, 
			     node1_t *PrevUMinI, node1_t *PrevVMinJ,
			     node1_t *UHead, process_variables *glob);
#if DEBUG_LEVEL > 0
void printSolution(process_variables *glob);
#endif


/******************************************************************************
float emd(signature_t *Signature1, signature_t *Signature2,
	  float (*Dist)(feature_t *, feature_t *), flow_t *Flow, int *FlowSize)
  
where

   Signature1, Signature2  Pointers to signatures that their distance we want
              to compute.
   Dist       Pointer to the ground distance. i.e. the function that computes
              the distance between two features.
   Flow       (Optional) Pointer to a vector of flow_t (defined in emd.h) 
              where the resulting flow will be stored. Flow must have n1+n2-1
              elements, where n1 and n2 are the sizes of the two signatures
              respectively.
              If NULL, the flow is not returned.
   FlowSize   (Optional) Pointer to an integer where the number of elements in
              Flow will be stored
              
******************************************************************************/

float emd(int size, double hist1[], double hist2[], double Dist[], double max_dist,
	  flow_t *Flow, int *FlowSize)
{
  int itr;
  double totalCost;
  double w;
  node2_t *XP;
  flow_t *FlowP;
  node1_t *U, *V;
  process_variables *glob;

  glob=(process_variables *)malloc(sizeof(process_variables));
  U=(node1_t*)malloc(MAX_SIG_SIZE1*sizeof(node1_t));
  V=(node1_t*)malloc(MAX_SIG_SIZE1*sizeof(node1_t));

  w = init(size, hist1, hist2, Dist, max_dist, glob);

#if DEBUG_LEVEL > 1
  printf("\nINITIAL SOLUTION:\n");
  printSolution(glob);
#endif
 
  if (glob->_n1 > 1 && glob->_n2 > 1)  /* IF _n1 = 1 OR _n2 = 1 THEN WE ARE DONE */
    {
      for (itr = 1; itr < MAX_ITERATIONS; itr++)
	{
	  /* FIND BASIC VARIABLES */
	  findBasicVariables(U, V, glob);
	  
	  /* CHECK FOR OPTIMALITY */
	  if (isOptimal(U, V,glob))
	    break;
	  
	  /* IMPROVE SOLUTION */
	  newSol(glob);
	  
#if DEBUG_LEVEL > 1
	  printf("\nITERATION # %d \n", itr);
	  printSolution(glob);
#endif
	}

//      if (itr == MAX_ITERATIONS)
//	G_message(_("emd: Maximum number of iterations has been reached (%d)\n"),MAX_ITERATIONS);
    }

  /* COMPUTE THE TOTAL FLOW */
  totalCost = 0;
  if (Flow != NULL)
    FlowP = Flow;
  for(XP=glob->_X; XP < glob->_EndX; XP++)
    {
      if (XP == glob->_EnterX)  /* _EnterX IS THE EMPTY SLOT */
	continue;
      if (XP->i == glob->_n1 || XP->j == glob->_n2)  /* DUMMY FEATURE */
	continue;
      if (XP->val == 0)  /* ZERO FLOW */
	continue;

      totalCost += (double)XP->val * glob->_C[NCINDEX(XP->i,XP->j,glob->_n1)];
      if (Flow != NULL)
	{
	  FlowP->from = XP->i;
	  FlowP->to = XP->j;
	  FlowP->amount = XP->val;
	  FlowP++;
	}
    }
  if (Flow != NULL)
    *FlowSize = FlowP-Flow;

#if DEBUG_LEVEL > 0
  printf("\n*** OPTIMAL SOLUTION (%d ITERATIONS): %f ***\n", itr, totalCost);
#endif

  free(U);
  free(V);
  free(glob);

  /* RETURN THE NORMALIZED COST == EMD */
  return (float)(totalCost / w);
}


/**********************
   init
**********************/
double init(int size, double hist1[], double hist2[], double Dist[], double max_dist, process_variables *glob)
{
  int i, j;
  double sSum, dSum, diff, *p;
  double *S, *D;

  S=(double *)malloc(MAX_SIG_SIZE1*sizeof(double));
  D=(double *)malloc(MAX_SIG_SIZE1*sizeof(double));
 
  glob->_n1 = size;
  glob->_n2 = size;

  if (glob->_n1 > MAX_SIG_SIZE)
      G_fatal_error(_("emd: Signature size is limited to %d\n"), MAX_SIG_SIZE);
  
  /* SET THE DISTANCE MATRIX */
  glob->_C = Dist;
  glob->_maxC = max_dist;
	
  /* SUM UP THE SUPPLY AND DEMAND */
  sSum = 0.0;
  for(i=0, p=hist1; i < glob->_n1; i++, p++)
    {
      S[i] = *p;
      sSum += *p;
      glob->_RowsX[i] = NULL;
    }
  dSum = 0.0;
  for(j=0, p=hist2; j < glob->_n2; j++, p++)
    {
      D[j] = *p;
      dSum += *p;
      glob->_ColsX[j] = NULL;
    }

  /* IF SUPPLY DIFFERENT THAN THE DEMAND, ADD A ZERO-COST DUMMY CLUSTER */
  diff = sSum - dSum;
  if (fabs(diff) >= EPSILON * sSum)
    {
      if (diff < 0.0)
	{
	  for (j=0; j < glob->_n2; j++)
	    glob->_C[NCINDEX(glob->_n1,j,glob->_n1)] = 0;
	  S[glob->_n1] = -diff;
	  glob->_RowsX[glob->_n1] = NULL;
	  glob->_n1++;
	}
      else
	{
	  for (i=0; i < glob->_n1; i++)
	    glob->_C[NCINDEX(i,glob->_n2,glob->_n1)] = 0;
	  D[glob->_n2] = diff;
	  glob->_ColsX[glob->_n2] = NULL;
	  glob->_n2++;
	}
    }

  /* INITIALIZE THE BASIC VARIABLE STRUCTURES */
  for (i=0; i < glob->_n1; i++)
    for (j=0; j < glob->_n2; j++)
	glob->_IsX[i][j] = 0;
  glob->_EndX = glob->_X;
   
  glob->_maxW = sSum > dSum ? sSum : dSum;

  /* FIND INITIAL SOLUTION */
  russel(S, D, glob);

  glob->_EnterX = glob->_EndX++;  /* AN EMPTY SLOT (ONLY _n1+_n2-1 BASIC VARIABLES) */

  free(S);
  free(D);

  return sSum > dSum ? dSum : sSum;
}


/**********************
    findBasicVariables
 **********************/
void findBasicVariables(node1_t *U, node1_t *V, process_variables *glob)
{
  int i, j, found;
  int UfoundNum, VfoundNum;
  node1_t u0Head, u1Head, *CurU, *PrevU;
  node1_t v0Head, v1Head, *CurV, *PrevV;

  /* INITIALIZE THE ROWS LIST (U) AND THE COLUMNS LIST (V) */
  u0Head.Next = CurU = U;
  for (i=0; i < glob->_n1; i++)
    {
      CurU->i = i;
      CurU->Next = CurU+1;
      CurU++;
    }
  (--CurU)->Next = NULL;
  u1Head.Next = NULL;

  CurV = V+1;
  v0Head.Next = glob->_n2 > 1 ? V+1 : NULL;
  for (j=1; j < glob->_n2; j++)
    {
      CurV->i = j;
      CurV->Next = CurV+1;
      CurV++;
    }
  (--CurV)->Next = NULL;
  v1Head.Next = NULL;

  /* THERE ARE _n1+_n2 VARIABLES BUT ONLY _n1+_n2-1 INDEPENDENT EQUATIONS,
     SO SET V[0]=0 */
  V[0].i = 0;
  V[0].val = 0;
  v1Head.Next = V;
  v1Head.Next->Next = NULL;

  /* LOOP UNTIL ALL VARIABLES ARE FOUND */
  UfoundNum=VfoundNum=0;
  while (UfoundNum < glob->_n1 || VfoundNum < glob->_n2)
    {

#if DEBUG_LEVEL > 3
      printf("UfoundNum=%d/%d,VfoundNum=%d/%d\n",UfoundNum,glob->_n1,VfoundNum,glob->_n2);
      printf("U0=");
      for(CurU = u0Head.Next; CurU != NULL; CurU = CurU->Next)
	printf("[%d]",CurU-U);
      printf("\n");
      printf("U1=");
      for(CurU = u1Head.Next; CurU != NULL; CurU = CurU->Next)
	printf("[%d]",CurU-U);
      printf("\n");
      printf("V0=");
      for(CurV = v0Head.Next; CurV != NULL; CurV = CurV->Next)
	printf("[%d]",CurV-V);
      printf("\n");
      printf("V1=");
      for(CurV = v1Head.Next; CurV != NULL; CurV = CurV->Next)
	printf("[%d]",CurV-V);
      printf("\n\n");
#endif
      
      found = 0;
      if (VfoundNum < glob->_n2)
	{
	  /* LOOP OVER ALL MARKED COLUMNS */
	  PrevV = &v1Head;
	  for (CurV=v1Head.Next; CurV != NULL; CurV=CurV->Next)
	    {
	      j = CurV->i;
	      /* FIND THE VARIABLES IN COLUMN j */
	      PrevU = &u0Head;
	      for (CurU=u0Head.Next; CurU != NULL; CurU=CurU->Next)
		{
		  i = CurU->i;
		  if (glob->_IsX[i][j])
		    {
		      /* COMPUTE U[i] */
		      CurU->val = glob->_C[NCINDEX(i,j,glob->_n1)] - CurV->val;
		      /* ...AND ADD IT TO THE MARKED LIST */
		      PrevU->Next = CurU->Next;
		      CurU->Next = u1Head.Next != NULL ? u1Head.Next : NULL;
		      u1Head.Next = CurU;
		      CurU = PrevU;
		    }
		  else
		    PrevU = CurU;
		}
	      PrevV->Next = CurV->Next;
	      VfoundNum++;
	      found = 1;
	    }
	}
     if (UfoundNum < glob->_n1)
	{
	  /* LOOP OVER ALL MARKED ROWS */
	  PrevU = &u1Head;
	  for (CurU=u1Head.Next; CurU != NULL; CurU=CurU->Next)
	    {
	      i = CurU->i;
	      /* FIND THE VARIABLES IN ROWS i */
	      PrevV = &v0Head;
	      for (CurV=v0Head.Next; CurV != NULL; CurV=CurV->Next)
		{
		  j = CurV->i;
		  if (glob->_IsX[i][j])
		    {
		      /* COMPUTE V[j] */
		      CurV->val = glob->_C[NCINDEX(i,j,glob->_n1)] - CurU->val;
		      /* ...AND ADD IT TO THE MARKED LIST */
		      PrevV->Next = CurV->Next;
		      CurV->Next = v1Head.Next != NULL ? v1Head.Next: NULL;
		      v1Head.Next = CurV;
		      CurV = PrevV;
		    }
		  else
		    PrevV = CurV;
		}
	      PrevU->Next = CurU->Next;
	      UfoundNum++;
	      found = 1;
	    }
	}
    if (! found)
	 G_fatal_error(_("emd: Unexpected error in findBasicVariables!\nThis typically happens when the EPSILON defined in\nemd.h is not right for the scale of the problem.\n"));
    }
}




/**********************
    isOptimal
 **********************/
int isOptimal(node1_t *U, node1_t *V, process_variables *glob)
{    
  double delta, deltaMin;
  int i, j, minI, minJ;

  /* FIND THE MINIMAL Cij-Ui-Vj OVER ALL i,j */
  deltaMin = INFINITY;
  for(i=0; i < glob->_n1; i++)
    for(j=0; j < glob->_n2; j++)
      if (! glob->_IsX[i][j])
	{
	  delta = glob->_C[NCINDEX(i,j,glob->_n1)] - U[i].val - V[j].val;
	  if (deltaMin > delta)
	    {
              deltaMin = delta;
	      minI = i;
	      minJ = j;
	    }
	}

#if DEBUG_LEVEL > 3
  printf("deltaMin=%f\n", deltaMin);
#endif

   if (deltaMin == INFINITY)
       G_fatal_error(_("emd: Unexpected error in isOptimal.\n"));

   
   glob->_EnterX->i = minI;
   glob->_EnterX->j = minJ;
   
   /* IF NO NEGATIVE deltaMin, WE FOUND THE OPTIMAL SOLUTION */
   return deltaMin >= -EPSILON * glob->_maxC;

/*
   return deltaMin >= -EPSILON;
 */
}


/**********************
    newSol
**********************/
void newSol(process_variables *glob)
{
    int i, j, k;
    double xMin;
    int steps;
    node2_t **Loop, *CurX, *LeaveX;

    Loop = (node2_t **)malloc(2*MAX_SIG_SIZE1*sizeof(node2_t*));
 
#if DEBUG_LEVEL > 3
    printf("EnterX = (%d,%d)\n", glob->_EnterX->i, glob->_EnterX->j);
#endif

    /* ENTER THE NEW BASIC VARIABLE */
    i = glob->_EnterX->i;
    j = glob->_EnterX->j;
    glob->_IsX[i][j] = 1;
    glob->_EnterX->NextC = glob->_RowsX[i];
    glob->_EnterX->NextR = glob->_ColsX[j];
    glob->_EnterX->val = 0;
    glob->_RowsX[i] = glob->_EnterX;
    glob->_ColsX[j] = glob->_EnterX;

    /* FIND A CHAIN REACTION */
    steps = findLoop(Loop, glob);

    /* FIND THE LARGEST VALUE IN THE LOOP */
    xMin = INFINITY;
    for (k=1; k < steps; k+=2)
      {
	if (Loop[k]->val < xMin)
	  {
	    LeaveX = Loop[k];
	    xMin = Loop[k]->val;
	  }
      }

    /* UPDATE THE LOOP */
    for (k=0; k < steps; k+=2)
      {
	Loop[k]->val += xMin;
	Loop[k+1]->val -= xMin;
      }

#if DEBUG_LEVEL > 3
    printf("LeaveX = (%d,%d)\n", LeaveX->i, LeaveX->j);
#endif

    /* REMOVE THE LEAVING BASIC VARIABLE */
    i = LeaveX->i;
    j = LeaveX->j;
    glob->_IsX[i][j] = 0;
    if (glob->_RowsX[i] == LeaveX)
      glob->_RowsX[i] = LeaveX->NextC;
    else
      for (CurX=glob->_RowsX[i]; CurX != NULL; CurX = CurX->NextC)
	if (CurX->NextC == LeaveX)
	  {
	    CurX->NextC = CurX->NextC->NextC;
	    break;
	  }
    if (glob->_ColsX[j] == LeaveX)
      glob->_ColsX[j] = LeaveX->NextR;
    else
      for (CurX=glob->_ColsX[j]; CurX != NULL; CurX = CurX->NextR)
	if (CurX->NextR == LeaveX)
	  {
	    CurX->NextR = CurX->NextR->NextR;
	    break;
	  }
    /* SET _EnterX TO BE THE NEW EMPTY SLOT */
    glob->_EnterX = LeaveX;

    free(Loop);
}



/**********************
    findLoop
**********************/
int findLoop(node2_t **Loop, process_variables *glob)
{
  int i, steps;
  node2_t **CurX, *NewX;
  char *IsUsed; 

  IsUsed=(char *)malloc(2*MAX_SIG_SIZE1*sizeof(char));
 
  for (i=0; i < glob->_n1+glob->_n2; i++)
    IsUsed[i] = 0;

  CurX = Loop;
  NewX = *CurX = glob->_EnterX;
  IsUsed[glob->_EnterX-glob->_X] = 1;
  steps = 1;

  do
    {
      if (steps%2 == 1)
	{
	  /* FIND AN UNUSED X IN THE ROW */
	  NewX = glob->_RowsX[NewX->i];
	  while (NewX != NULL && IsUsed[NewX-glob->_X])
	    NewX = NewX->NextC;
	}
      else
	{
	  /* FIND AN UNUSED X IN THE COLUMN, OR THE ENTERING X */
	  NewX = glob->_ColsX[NewX->j];
	  while (NewX != NULL && IsUsed[NewX-glob->_X] && NewX != glob->_EnterX)
	    NewX = NewX->NextR;
	  if (NewX == glob->_EnterX)
	    break;
 	}

     if (NewX != NULL)  /* FOUND THE NEXT X */
       {
	 /* ADD X TO THE LOOP */
	 *++CurX = NewX;
	 IsUsed[NewX-glob->_X] = 1;
	 steps++;
#if DEBUG_LEVEL > 3
	 printf("steps=%d, NewX=(%d,%d)\n", steps, NewX->i, NewX->j);    
#endif
       }
     else  /* DIDN'T FIND THE NEXT X */
       {
	 /* BACKTRACK */
	 do
	   {
	     NewX = *CurX;
	     do 
	       {
		 if (steps%2 == 1)
		   NewX = NewX->NextR;
		 else
		   NewX = NewX->NextC;
	       } while (NewX != NULL && IsUsed[NewX-glob->_X]);
	     
	     if (NewX == NULL)
	       {
		 IsUsed[*CurX-glob->_X] = 0;
		 CurX--;
		 steps--;
	       }
	 } while (NewX == NULL && CurX >= Loop);
	 
#if DEBUG_LEVEL > 3
	 printf("BACKTRACKING TO: steps=%d, NewX=(%d,%d)\n",
		steps, NewX->i, NewX->j);    
#endif
           IsUsed[*CurX-glob->_X] = 0;
	   *CurX = NewX;
	   IsUsed[NewX-glob->_X] = 1;
       }     
    } while(CurX >= Loop);
  
  if (CurX == Loop)
      G_fatal_error(_("emd: Unexpected error in findLoop!\n"));

#if DEBUG_LEVEL > 3
  printf("FOUND LOOP:\n");
  for (i=0; i < steps; i++)
    printf("%d: (%d,%d)\n", i, Loop[i]->i, Loop[i]->j);
#endif

  free(IsUsed);
  return steps;
}



/**********************
    russel
**********************/
void russel(double *S, double *D, process_variables *glob)
{
  int i, j, found, minI, minJ;
  double deltaMin, oldVal, diff;
  double *Delta;
  node1_t *Ur, *Vr;
  node1_t uHead, *CurU, *PrevU;
  node1_t vHead, *CurV, *PrevV;
  node1_t *PrevUMinI, *PrevVMinJ, *Remember;

  Delta=(double*)malloc(MAX_SIG_SIZE1*MAX_SIG_SIZE1*sizeof(double));
  Ur=(node1_t*)malloc(MAX_SIG_SIZE1*sizeof(node1_t));
  Vr=(node1_t*)malloc(MAX_SIG_SIZE1*sizeof(node1_t));

  /* INITIALIZE THE ROWS LIST (Ur), AND THE COLUMNS LIST (Vr) */
  uHead.Next = CurU = Ur;
  for (i=0; i < glob->_n1; i++)
    {
      CurU->i = i;
      CurU->val = -INFINITY;
      CurU->Next = CurU+1;
      CurU++;
    }
  (--CurU)->Next = NULL;
  
  vHead.Next = CurV = Vr;
  for (j=0; j < glob->_n2; j++)
    {
      CurV->i = j;
      CurV->val = -INFINITY;
      CurV->Next = CurV+1;
      CurV++;
    }
  (--CurV)->Next = NULL;
  
  /* FIND THE MAXIMUM ROW AND COLUMN VALUES (Ur[i] AND Vr[j]) */
  for(i=0; i < glob->_n1 ; i++)
    for(j=0; j < glob->_n2 ; j++)
      {
	float v;
	v = glob->_C[NCINDEX(i,j,glob->_n1)];
	if (Ur[i].val <= v)
	  Ur[i].val = v;
	if (Vr[j].val <= v)
	  Vr[j].val = v;
      }
  
  /* COMPUTE THE Delta MATRIX */
  for(i=0; i < glob->_n1 ; i++)
    for(j=0; j < glob->_n2 ; j++)
      Delta[NCINDEX(i,j,MAX_SIG_SIZE1)] = glob->_C[NCINDEX(i,j,glob->_n1)] - Ur[i].val - Vr[j].val;

  /* FIND THE BASIC VARIABLES */
  do
    {
#if DEBUG_LEVEL > 3
      printf("Ur=");
      for(CurU = uHead.Next; CurU != NULL; CurU = CurU->Next)
	printf("[%d]",CurU-Ur);
      printf("\n");
      printf("Vr=");
      for(CurV = vHead.Next; CurV != NULL; CurV = CurV->Next)
	printf("[%d]",CurV-Vr);
      printf("\n");
      printf("\n\n");
#endif
 
      /* FIND THE SMALLEST Delta[i][j] */
      found = 0; 
      deltaMin = INFINITY;      
      PrevU = &uHead;
      for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
	{
	  int i;
	  i = CurU->i;
	  PrevV = &vHead;
	  for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
	    {
	      int j;
	      j = CurV->i;
	      if (deltaMin > Delta[NCINDEX(i,j,MAX_SIG_SIZE1)])
		{
		  deltaMin = Delta[NCINDEX(i,j,MAX_SIG_SIZE1)];
		  minI = i;
		  minJ = j;
		  PrevUMinI = PrevU;
		  PrevVMinJ = PrevV;
		  found = 1;
		}
	      PrevV = CurV;
	    }
	  PrevU = CurU;
	}
      
      if (! found)
	break;

      /* ADD X[minI][minJ] TO THE BASIS, AND ADJUST SUPPLIES AND COST */
      Remember = PrevUMinI->Next;
      addBasicVariable(minI, minJ, S, D, PrevUMinI, PrevVMinJ, &uHead, glob);

      /* UPDATE THE NECESSARY Delta[][] */
      if (Remember == PrevUMinI->Next)  /* LINE minI WAS DELETED */
	{
	  for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
	    {
	      int j;
	      j = CurV->i;
	      if (CurV->val == glob->_C[NCINDEX(minI,j,glob->_n1)])  /* COLUMN j NEEDS UPDATING */
		{
		  /* FIND THE NEW MAXIMUM VALUE IN THE COLUMN */
		  oldVal = CurV->val;
		  CurV->val = -INFINITY;
		  for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
		    {
		      int i;
		      i = CurU->i;
		      if (CurV->val <= glob->_C[NCINDEX(i,j,glob->_n1)])
			CurV->val = glob->_C[NCINDEX(i,j,glob->_n1)];
		    }
		  
		  /* IF NEEDED, ADJUST THE RELEVANT Delta[*][j] */
		  diff = oldVal - CurV->val;
		  if (fabs(diff) < EPSILON * glob->_maxC)
		    for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
		      Delta[NCINDEX(CurU->i,j,MAX_SIG_SIZE1)] += diff;
		}
	    }
	}
      else  /* COLUMN minJ WAS DELETED */
	{
	  for (CurU=uHead.Next; CurU != NULL; CurU=CurU->Next)
	    {
	      int i;
	      i = CurU->i;
	      if (CurU->val == glob->_C[NCINDEX(i,minJ,glob->_n1)])  /* ROW i NEEDS UPDATING */
		{
		  /* FIND THE NEW MAXIMUM VALUE IN THE ROW */
		  oldVal = CurU->val;
		  CurU->val = -INFINITY;
		  for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
		    {
		      int j;
		      j = CurV->i;
		      if(CurU->val <= glob->_C[NCINDEX(i,j,glob->_n1)])
			CurU->val = glob->_C[NCINDEX(i,j,glob->_n1)];
		    }
		  
		  /* If NEEDED, ADJUST THE RELEVANT Delta[i][*] */
		  diff = oldVal - CurU->val;
		  if (fabs(diff) < EPSILON * glob->_maxC)
		    for (CurV=vHead.Next; CurV != NULL; CurV=CurV->Next)
		      Delta[NCINDEX(i,CurU->i,MAX_SIG_SIZE1)] += diff;
		}
	    }
	}
    } while (uHead.Next != NULL || vHead.Next != NULL);
    free(Delta);
    free(Ur);
    free(Vr);
}




/**********************
    addBasicVariable
**********************/
void addBasicVariable(int minI, int minJ, double *S, double *D, 
			     node1_t *PrevUMinI, node1_t *PrevVMinJ,
			     node1_t *UHead, process_variables *glob)
{
  double T;
  
  if (fabs(S[minI]-D[minJ]) <= EPSILON * glob->_maxW)  /* DEGENERATE CASE */
    {
      T = S[minI];
      S[minI] = 0;
      D[minJ] -= T; 
    }
  else if (S[minI] < D[minJ])  /* SUPPLY EXHAUSTED */
    {
      T = S[minI];
      S[minI] = 0;
      D[minJ] -= T; 
    }
  else  /* DEMAND EXHAUSTED */
    {
      T = D[minJ];
      D[minJ] = 0; 
      S[minI] -= T; 
    }

  /* X(minI,minJ) IS A BASIC VARIABLE */
  glob->_IsX[minI][minJ] = 1; 

  glob->_EndX->val = T;
  glob->_EndX->i = minI;
  glob->_EndX->j = minJ;
  glob->_EndX->NextC = glob->_RowsX[minI];
  glob->_EndX->NextR = glob->_ColsX[minJ];
  glob->_RowsX[minI] = glob->_EndX;
  glob->_ColsX[minJ] = glob->_EndX;
  glob->_EndX++;

  /* DELETE SUPPLY ROW ONLY IF THE EMPTY, AND IF NOT LAST ROW */
  if (S[minI] == 0 && UHead->Next->Next != NULL)
    PrevUMinI->Next = PrevUMinI->Next->Next;  /* REMOVE ROW FROM LIST */
  else
    PrevVMinJ->Next = PrevVMinJ->Next->Next;  /* REMOVE COLUMN FROM LIST */
}





/**********************
    printSolution
**********************/
void printSolution(process_variables *glob)
{
  node2_t *P;
  double totalCost;

  totalCost = 0;

#if DEBUG_LEVEL > 2
  printf("SIG1\tSIG2\tFLOW\tCOST\n");
#endif
  for(P=glob->_X; P < glob->_EndX; P++)
    if (P != glob->_EnterX && glob->_IsX[P->i][P->j])
      {
#if DEBUG_LEVEL > 2
	printf("%d\t%d\t%f\t%f\n", P->i, P->j, P->val, glob->_C[NCINDEX(P->i,P->j,glob->_n1)]);
#endif
	totalCost += (double)P->val * glob->_C[NCINDEX(P->i,P->j,glob->_n1)];
      }

  printf("COST = %f\n", totalCost);
}


