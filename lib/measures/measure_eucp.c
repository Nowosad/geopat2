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
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>


double **euclidean_period_create_mtx(double *s1, double *s2, int dim, int *dims) {
  int i,j,k,n,d;
  double sum, p;
  double **mtx;
  double *h1,*h2;

  if(dim==1) {
    n = *dims;
    d = 1;
  } else {
    n = dims[1];
    d = dims[0];
  }
  mtx = (double **)calloc(n, sizeof(double *));

  for(i=0; i<n; i++)
    mtx[i] = (double *)malloc(n*sizeof(double));

  for(i=0; i<n; i++)
    for(j=i; j<n; j++) {
      sum = 0.0;
      h1 = s1 + i;
      h2 = s2 + j;
      for(k=0; k<d; k++) {
        p = (*h1) - (*h2);
        sum += p*p;
        h1 += n;
        h2 += n;
      }
      mtx[i][j] = sum;
    }

  return mtx;
}

void euclidean_period_free_mtx(double **mtx, int n) {
  int i;
  for(i=0; i<n; i++)
    free(mtx[i]);
  free(mtx);
}

double euclidean_period_find_min_sum(double **mtx, int n){
  int i,j,s;
  double sum, min_sum = DBL_MAX;

  for(s=0; s<n; s++) {
    sum = 0.0;
    for(i=0; i<n; i++) {
      j = (i+s)%n;
      sum += (i<j)?mtx[i][j]:mtx[j][i];
    }
    if(sum<min_sum)
      min_sum = sum;
  }
  return min_sum;
}


double euclidean_period(double **signatures, int num_of_signatures, int size_of_signature, int num_dims, int* dims, ...) {

  double distance;
  double **mtx;

  if(num_of_signatures<2 || num_dims>2) return 0.0;

  mtx = euclidean_period_create_mtx(signatures[0],signatures[1],num_dims,dims);
  if(num_dims==1) {
    distance = euclidean_period_find_min_sum(mtx,dims[0]);
    euclidean_period_free_mtx(mtx,dims[0]);
  } else {
    distance = euclidean_period_find_min_sum(mtx,dims[1]);
    euclidean_period_free_mtx(mtx,dims[1]);
  }

  return sqrt(distance);
}

