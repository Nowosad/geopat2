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
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>

#include <stdio.h>

/*
#define MAXIMUM2(a,b) ((a) > (b) ? (a):(b))
#define MINIMUM2(a,b) ((a) < (b) ? (a):(b))
#define MINIMUM3(a,b,c) (MINIMUM2(MINIMUM2((a),(b)),(c)))
*/

double *_DTW_create_mtx1(int n) {
	double *p = malloc(n*n*sizeof(double));
	assert(p!=NULL);
	return p;
}

double **_DTW_create_mtx2(double *d, int n) {\
	int i;
	double **p = malloc(n*sizeof(double*));
	assert(p!=NULL);
	for(i=0; i<n; i++) {
		p[i]=d+n*i;
	}
	return p;
}

void _DTW_calc_data_distance2(int n, int dim, double *v1, double *v2, double **d) {
	int i, j, k;
	double *p1, *p2;
	double pd, dist;

	for(i=0; i<n; i++)
	  for(j=0; j<n; j++) {
            p1=v1+i;
	    p2=v2+j;
	    pd=*p1-*p2;
            dist=pd*pd;
		for(k=1; k<dim; k++) {
		  p1+=n;
		  p2+=n;
		  pd=*p1-*p2;
		  dist+=pd*pd;
		}
		d[i][j]=dist;
	  }
}

void _DTW_rotate_data_distances(double **d, int n) {
	int i;
	double *p;
	
	p=d[0];
	for(i=1; i<n; i++) {
		d[i-1]=d[i];
    }
    d[n-1]=p;
}

double _DTW_dtw_min(double x,double y,double z) {
    if((x<=y) && (x<=z)) return x;
    if((y<=x) && (y<=z)) return y;
    if((z<=x) && (z<=y)) return z;
    return 0.0;
}

double _DTW_dtw_distance(double **cd, double **d, int n) {
	int i, j;
	
    cd[0][0]=d[0][0];
    for(i=1; i<n; i++) {
      cd[i][0]=cd[i-1][0] + d[i][0];
      cd[0][i]=cd[0][i-1] + d[0][i];
    }

    for(i=1; i<n; i++)
      for(j=1; j<n; j++)
        cd[i][j]=d[i][j]+_DTW_dtw_min(cd[i-1][j],cd[i][j-1],cd[i-1][j-1]);

    return cd[n-1][n-1];
}

double _DTW_euc_distance(double **d, int n) {
	double dist;
	int i;

    dist=0.0;
    for(i=0; i<n; i++) 
        dist+=d[i][i];

    return dist;
}

int _DTW_dtw_trace(double **cd, int n) {
  int i, j, len;

    i = n-1;
    j = n-1;
//printf("======== %d =========\n",n);
    len = 1;
    while ((i > 0) || (j > 0)) {
//printf(" step %d - (%d, %d)\n",len,i,j);
      if(i == 0)
        j--;
      else if(j == 0)
        i--;
      else if(cd[i-1][j] == _DTW_dtw_min(cd[i-1][j],cd[i][j-1],cd[i-1][j-1]))
        i--;
      else if(cd[i][j-1] == _DTW_dtw_min(cd[i-1][j],cd[i][j-1],cd[i-1][j-1]))
        j--;
      else {
        i--;
        j--;
      }
      len++;
    }
//printf(" step %d - (%d, %d)\n",len,i,j);
  return len;
}


double tsDTW(double **signatures, int num_of_signatures, int size_of_signature, int num_of_dims, int* dims, ...) { 

  double *m_d;
  double **d;
  double *m_cd;
  double **cd;
  double dist,len,scale=1.0;
  int n, dim;

  if(num_of_dims<1 || num_of_dims>2) return 0.0;
  if(num_of_signatures<2) return 0.0;

  n=dims[0];
  if(num_of_dims==1) dim = 1; else  { dim = dims[0]; n=dims[1];}

  m_d = _DTW_create_mtx1(n);
  d = _DTW_create_mtx2(m_d,n);
  m_cd = _DTW_create_mtx1(n);
  cd = _DTW_create_mtx2(m_cd,n);

  _DTW_calc_data_distance2(n,dim,signatures[0],signatures[1],d);
  dist=_DTW_dtw_distance(cd,d,n);
  len=_DTW_dtw_trace(cd,n);
  free(cd);
  free(d);
  free(m_cd);
  free(m_d);
  

  if(dim>1)
    scale=1.0/sqrt(dim);
  return scale*sqrt(dist)/len;
//  return scale*sqrt(dist);
//  return sqrt(dist);
}

double tsDTWP(double **signatures, int num_of_signatures, int size_of_signature, int num_of_dims, int* dims, ...) { 

  double *m_d;
  double **d;
  double *m_cd;
  double **cd;
  double dist,pd,scale=1.0;
  int n, dim, i, len;

  if(num_of_dims<1 || num_of_dims>2) return 0.0;
  if(num_of_signatures<2) return 0.0;

  n=dims[0];
  if(num_of_dims==1) dim = 1; else  { dim = dims[0]; n=dims[1];}

  m_d = _DTW_create_mtx1(n);
  d = _DTW_create_mtx2(m_d,n);
  m_cd = _DTW_create_mtx1(n);
  cd = _DTW_create_mtx2(m_cd,n);

  _DTW_calc_data_distance2(n,dim,signatures[0],signatures[1],d);
  dist=_DTW_dtw_distance(cd,d,n);
  len=_DTW_dtw_trace(cd,n);
  for(i=1; i<n; i++) {
    pd=_DTW_dtw_distance(cd,d,n);
    if(pd<dist) {
	  dist=pd;
      len=_DTW_dtw_trace(cd,n);
    }
    _DTW_rotate_data_distances(d,n); 
  }
  free(cd);
  free(d);
  free(m_cd);
  free(m_d);

  if(dim>1)
    scale=1.0/sqrt(dim);
//  return scale*sqrt(dist)/len;
  return scale*sqrt(dist);
//  return sqrt(dist);
}


double tsEUC(double **signatures, int num_of_signatures, int size_of_signature, int num_of_dims, int* dims, ...) { 

  double *m_d;
  double **d;
  double dist,scale=1.0;
  int n, dim;

  if(num_of_dims<1 || num_of_dims>2) return 0.0;
  if(num_of_signatures<2) return 0.0;

  n=dims[0];
  if(num_of_dims==1) dim = 1; else { dim = dims[0]; n=dims[1];}

  m_d = _DTW_create_mtx1(n);
  d = _DTW_create_mtx2(m_d,n);

  _DTW_calc_data_distance2(n,dim,signatures[0],signatures[1],d);
  dist=_DTW_euc_distance(d,n);
  free(d);
  free(m_d);

  if(dim>1)
    scale=1.0/sqrt(dim);
  return scale*sqrt(dist);
}

double tsEUCP(double **signatures, int num_of_signatures, int size_of_signature, int num_of_dims, int* dims, ...) { 

  double *m_d;
  double **d;
  double dist,pd,scale=1.0;
  int n, dim, i;

  if(num_of_dims<1 || num_of_dims>2) return 0.0;
  if(num_of_signatures<2) return 0.0;

  n=dims[0];
  if(num_of_dims==1) dim = 1; else { dim = dims[0]; n=dims[1];}

  m_d = _DTW_create_mtx1(n);
  d = _DTW_create_mtx2(m_d,n);

  _DTW_calc_data_distance2(n,dim,signatures[0],signatures[1],d);
  dist=_DTW_euc_distance(d,n);
  for(i=1; i<n; i++) {
    pd=_DTW_euc_distance(d,n);
    if(pd<dist)
	  dist=pd;
    _DTW_rotate_data_distances(d,n); 
  }
  free(d);
  free(m_d);

  if(dim>1)
    scale=1.0/sqrt(dim);
  return scale*sqrt(dist);
//  return sqrt(dist);
}

double tsDTWPa(double **signatures, int num_of_signatures, int size_of_signature, int num_of_dims, int* dims, ...) { 

  double *m_d;
  double **d;
  double **d_align;
  double *m_cd;
  double **cd;
  double dist,pd,scale=1.0;
  int n, dim, i, len, size;

  if(num_of_dims<1 || num_of_dims>2) return 0.0;
  if(num_of_signatures<2) return 0.0;

  n=dims[0];
  if(num_of_dims==1) dim = 1; else dim = dims[1];

  m_d = _DTW_create_mtx1(n);
  d = _DTW_create_mtx2(m_d,n);
  size = n*sizeof(double*);
  d_align = malloc(size);
  m_cd = _DTW_create_mtx1(n);
  cd = _DTW_create_mtx2(m_cd,n);

  _DTW_calc_data_distance2(n,dim,signatures[0],signatures[1],d);
  dist=_DTW_euc_distance(d,n);
  len=_DTW_dtw_trace(cd,n);
  for(i=1; i<n; i++) {
    pd=_DTW_euc_distance(d,n);
    if(pd<dist) {
	  dist=pd;
	  memcpy(d_align,d,size);
	}
    _DTW_rotate_data_distances(d,n); 
  }

  dist = _DTW_dtw_distance(cd,d_align,n);
  len = _DTW_dtw_trace(cd,n);
  free(cd);
  free(d);
  free(d_align);
  free(m_cd);
  free(m_d);

  if(dim>1)
    scale=1.0/sqrt(dim);
  return (float)(scale*sqrt(dist));
//  return (float)(scale*sqrt(dist)/len);
//  return (float)(sqrt(dist));
}
