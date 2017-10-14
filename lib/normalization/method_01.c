/****************************************************************************
 *
 * MODULE:	Signature normalization library
 * AUTHOR(S):	Pawel Netzel
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <float.h>

int normalization_01(double *signature, int signature_len) {

  int i;
  double min, max, *p;

  min = DBL_MAX;
  max = DBL_MIN;
  p = signature;
  for(i=0; i<signature_len; i++) {
    if(*p>max) max = *p;
    else if(*p<min) min = *p;
    p++;
  }

  if(max<=min) return 1;

  max = 1.0/(max-min);
  p = signature;
  for(i=0; i<signature_len; i++) {
    *p = (*p - min) * max;
    p++;
  }

  return 0;
}

