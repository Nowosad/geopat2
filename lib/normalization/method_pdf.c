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
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

int normalization_pdf(double *signature, int signature_len) {

  assert(signature!=NULL);

  int i;
  double sum, *p;

  sum = 0.0;
  p = signature;
  for(i=0; i<signature_len; i++) {
    sum += (*p);
    p++;
  }

  if(sum<0.0) return 1;
  if(sum==0.0) return 0;

  sum = 1.0/sum;
  p = signature;
  for(i=0; i<signature_len; i++) {
    (*p) *= sum;
    p++;
//printf(".");
  }
//printf("\n"); fflush(stdout);

  return 0;
}

