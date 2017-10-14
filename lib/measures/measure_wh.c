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

double wave_hedges(double **signatures, int num_of_signatures, int size_of_signature, int num_dims, int* dims, ...) {
  int i;

  double h0, h1;
  double dist = 0.0;
  double d    = 0.0;
  int cnt_non_empty = 0;

  if(num_of_signatures<2) return 0.0;

  for(i=0; i<size_of_signature; i++) {
    h0 = signatures[0][i];
    h1 = signatures[1][i];
    d  = (h0<h1)?h1:h0;
    if(d > 0) {
      h0 = h0-h1;
      dist += ((h0<0.0)?-h0:h0)/d;
      cnt_non_empty++;
    }
  }

  dist /= (double)cnt_non_empty;
  return dist;
}

