/****************************************************************************
 *
 * MODULE:	segmentation quality of motifels grid
 * AUTHOR(S):	Jaroslaw Jasiewicz, Jacek Niesterowicz, Tomasz Stepinski
 * PURPOSE:	information retrival using categorical maps:
 *		calculates quality of segmentation
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include "local_proto.h"

int sort_desc (const void * a, const void * b)
{
   return ( *(int*)b - *(int*)a );
   /* sort descending */
 }

int sort_asc (const void * a, const void * b)
{
   /* sort ascending */
   return ( *(int*)a - *(int*)b );
}

int sort_double_desc (const void * a, const void * b)
{
   /* sort descending */
   return (int)( 1000000*(*(double*)b) - 1000000*(*(double*)a) );
}

int sort_double_asc (const void * a, const void * b)
{
   /* sort ascending */
   return (int)( 1000000*(*(double*)a) - 1000000*(*(double*)b) );
}

/* ==================================================================== */

double* use_histogram(DATAINFO* d, int index)
{
	return d->all_histograms[index];
	/* this function will be extended in the future to get data from SSD or from memory*/
}

int add_histograms(DATAINFO* d, double* o, double* h, int num_of_areas)
{
	double o_div=(num_of_areas-1.)/(double)num_of_areas;
	double h_div=1./(double)num_of_areas;
	int i;
	for(i=0;i<d->size_of_histogram;++i)
		o[i]=o[i]*o_div+h[i]*h_div;
	return 0;
}
