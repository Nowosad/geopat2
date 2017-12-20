/****************************************************************************
 *
 * MODULE:	Similarity measures library (triangular measure)
 * AUTHOR(S):	Jaroslaw Jasiewicz, Jakub Nowosad
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <math.h>
#include <stdarg.h>

double triangular(double **signatures, int num_of_signatures, int size_of_signature, int num_dims, int* dims, ...)
{
        /* distance */
        /* euclidian distance normalised by sum of all signatures */
        int i;
        double dist     = 0.0;
        double divident = 0.0;
        double divisor  = 0.0;
        
        for(i=0; i<size_of_signature; i++) {
                divident = (signatures[0][i] - signatures[1][i]) * (signatures[0][i] - signatures[1][i]);
                divisor = signatures[0][i] + signatures[1][i];
                dist += (divisor>0)?divident/divisor:0;
        }
        dist = fabs(dist/2.0);
        /* dist = sqrt(dist); */
        /* dist = 1-dist; */
        return dist;
}
