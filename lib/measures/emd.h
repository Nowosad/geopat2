#ifndef _EMD_H
#define _EMD_H
/*
    emd.h

    Changes 12/15/2016: 
	adaptation to GeoPAT 2.0 (GRASS free)
    Author: Pawel Netzel, Space Informatics Lab

    Changes 11/05/2014:
	cost matrix transfer from similarity library,
	error handling from GRASS,
	using dynamically allocated arrays.
    Author: Pawel Netzel, Space Informatics Lab

    Last update: 3/24/98

    An implementation of the Earth Movers Distance.
    Based of the solution for the Transportation problem as described in
    "Introduction to Mathematical Programming" by F. S. Hillier and
    G. J. Lieberman, McGraw-Hill, 1990.

    Copyright (C) 1998 Yossi Rubner
    Computer Science Department, Stanford University
    E-Mail: rubner@cs.stanford.edu   URL: http://vision.stanford.edu/~rubner
*/


/* DEFINITIONS */
#define MAX_SIG_SIZE   2000
#define MAX_ITERATIONS 1000
//#define INFINITY       1e20 /* thing to remove, already defined */
#define EPSILON        1e-6


typedef struct
{
  int from;             /* Feature number in signature 1 */
  int to;               /* Feature number in signature 2 */
  float amount;         /* Amount of flow from "from" to "to" */
} flow_t;



float emd(int size, double hist1[], double hist2[],
	  double Dist[], double max_dist,
	  flow_t *Flow, int *FlowSize);

#endif
