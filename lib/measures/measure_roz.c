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

double rozicka(double **signatures, int num_of_signatures, int size_of_signature, int num_dims, int* dims, ...) {
  int i;
  double mi=0.0;
  double ma=0.0;
  double h0, h1, similarity;

  if(num_of_signatures<2) return 0.0;

  for(i=0; i<size_of_signature; i++) {
    h0 = signatures[0][i];
    h1 = signatures[1][i];
    mi+=(h0>h1)?h1:h0;
    ma+=(h0<h1)?h1:h0;
  }

  similarity = mi/ma;

  return 1.0-similarity;
}
