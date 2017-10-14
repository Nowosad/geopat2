/****************************************************************************
 *
 * MODULE:	Similarity measures library
 * AUTHOR(S):	Pawel Netzel
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <math.h>
#include <stdarg.h>

double jaccard(double **signatures, int num_of_signatures, int size_of_signature, int num_dims, int* dims, ...) {

  int i;

  double d = 0.0, d0 = 0.0, d1 = 0.0;

  if(num_of_signatures<2) return 0.0;

  for(i=0; i<size_of_signature; i++) {
    d  += signatures[0][i]*signatures[1][i];
    d0 += signatures[0][i]*signatures[0][i];
    d1 += signatures[1][i]*signatures[1][i];
  }

  return 1.0 - d/(d0+d1-d);
}

