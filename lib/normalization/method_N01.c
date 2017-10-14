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
#include <math.h>


int normalization_N01(double *signature, int signature_len) {

  int i;
  double n, avg, std, sumx, sumxx, *p;

  sumx = 0;
  sumxx = 0;
  p = signature;
  for(i=0; i<signature_len; i++) {
    sumx += *p;
    sumxx += (*p) * (*p);
    p++;
  }

  avg = sumx/signature_len;
  if(signature_len<2) 
    std = 1.0;
  else {
    n=1.0/(signature_len-1);
    std = sqrt(sumxx*n-2*n*avg*sumx+n*avg*avg);
  }

  if(std ==0) std = 1.0; else std = 1.0/std;
  p = signature;
  for(i=0; i<signature_len; i++) {
    *p = (*p - avg) * std;
    p++;
  }

  return 0;

}

