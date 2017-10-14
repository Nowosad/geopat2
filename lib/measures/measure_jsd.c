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

double jensen_shannon(double **signatures, int num_of_signatures, int size_of_signature, int num_dims, int* dims, ...) {
  double w = 1.0/(double)num_of_signatures;
  double sumE, sumH, entH, h;
  int i,j;
  double lg, distance = 0.0;

  if(num_of_signatures<2) return 0.0;

  lg = 1.0/log(2);
  
  for(j=0; j<size_of_signature; j++) {
    sumH = 0.0;
    sumE = 0.0;
    for(i=0; i<num_of_signatures; i++) {
      h=signatures[i][j];
      if(h>0) {
        sumH+=h;
	sumE+=h*lg*log(h);
      }
    }
    sumE *= w;
    sumH *= w;

    entH = (sumH==0)?0.0:sumH*lg*log(sumH);
    distance +=sumE-entH;
  }
  if(num_of_signatures>2)
    distance /= lg*log((double)num_of_signatures);

  if(distance>1.0) distance = 1.0;
  
  return distance;
}

