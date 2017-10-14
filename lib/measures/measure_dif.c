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

double ardiff(double **signatures, int num_of_signatures, int size_of_signature, int num_dims, int* dims, ...) { 
  int i;
  double h0, h1, distance = 0.0;

  if(num_of_signatures<2) return 0.0;

  for(i=0; i<size_of_signature; i++) {
	h0 = signatures[0][i];
	h1 = signatures[1][i];
    distance += ((h0<h1)?h1:h0)-((h0>h1)?h1:h0);
  }

  return sqrt(0.5*distance);
}

