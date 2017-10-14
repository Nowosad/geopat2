/****************************************************************************
 *
 * MODULE:	Coocurence matrix
 * AUTHOR(S):	Pawel Netzel
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <math.h>
#include <ezgdal.h>
#include <assert.h>
#include <stdarg.h>

int jcov_triangular_index(int r, int c) {
  r++; c++;
  if(c<=r)
    return (r-1)*r/2+c-1;
  else
    return (c-1)*c/2+r-1;

//  return (n*c1)+c2-((c1*(c1+1))/2);
}


int jcov(EZGDAL_FRAME **frames, int num_of_frames, double *signature, int signature_len, ...) {
  int i, r, c, cat1, cat2, cols, rows;
  EZGDAL_LAYER *l;
  EZGDAL_FRAME *f;
  double v;
  double *sx, *sy, *ssx, *ssy;
  int *N, sumN;
  double wx, wy, x, y;
  int cx, cy, len;

  f = frames[0];
  l = f->owner.stripe->layer;
  cols = f->cols;
  rows = f->rows;
  len = signature_len/2;

  cx = 1+cols/2;
  cy = 1+rows/2;
  wx = 0.5/(double)cx;
  wy = 0.5/(double)cy;
  sx  = (double *)calloc(signature_len, sizeof(double));
  ssx = (double *)calloc(signature_len, sizeof(double));
  sy  = (double *)calloc(signature_len, sizeof(double));
  ssy = (double *)calloc(signature_len, sizeof(double));
  N = (int *)calloc(signature_len, sizeof(int));

  for(r=0; r<rows-1; r++) {
    for(c=0; c<cols-1; c++) {

      v = f->buffer[r][c];
      if(!ezgdal_is_null(l,v)) {

        cat1 = ezgdal_get_value_index(l, v);
        assert(cat1>=0);

        v = f->buffer[r+1][c];
        if(!ezgdal_is_null(l,v)) {


          cat2 = ezgdal_get_value_index(l, v);
          assert(cat2>=0);

          i = jcov_triangular_index(cat1,cat2);
          x = wx*(c-cx);
          y = wy*(r-cy+0.5);
          sx[i]  += x;
          ssx[i] += x*x;
          sy[i]  += y;
          ssy[i] += y*y;
          N[i]++;
        }

        v = f->buffer[r][c+1];
        if(!ezgdal_is_null(l,v)) {

          cat2 = ezgdal_get_value_index(l, v);
          assert(cat2>=0);
          i = jcov_triangular_index(cat1,cat2);
          x = wx*(c-cx+0.5);
          y = wy*(r-cy);
          sx[i]  += x;
          ssx[i] += x*x;
          sy[i]  += y;
          ssy[i] += y*y;
          N[i]++;
        }
      }
    }
  }

  for(r=0; r<rows-1; r++) {

    v = f->buffer[r][cols-1];
    if(!ezgdal_is_null(l,v)) {

      cat1 = ezgdal_get_value_index(l, v);
      assert(cat1>=0);

      v = f->buffer[r+1][cols-1];
      if(!ezgdal_is_null(l,v)) {

        cat2 = ezgdal_get_value_index(l, v);
        assert(cat2>=0);

        i = jcov_triangular_index(cat1,cat2);
        x = wx*(c-cx);
        y = wy*(r-cy+0.5);
        sx[i]  += x;
        ssx[i] += x*x;
        sy[i]  += y;
        ssy[i] += y*y;
        N[i]++;
      }
    }
  }

  for(c=0; c<cols-1; c++) {

    v = f->buffer[rows-1][c];
    if(!ezgdal_is_null(l,v)) {

      cat1 = ezgdal_get_value_index(l, v);
      assert(cat1>=0);

      v = f->buffer[rows-1][c+1];
      if(!ezgdal_is_null(l,v)) {

        cat2 = ezgdal_get_value_index(l, v);
        assert(cat2>=0);

        i = jcov_triangular_index(cat1,cat2);
        x = wx*(c-cx+0.5);
        y = wy*(r-cy);
        sx[i]  += x;
        ssx[i] += x*x;
        sy[i]  += y;
        ssy[i] += y*y;
        N[i]++;
      }
    }
  }

  sumN = 0;
  for(i=0; i<len; i++)
    sumN += N[i];

  wy = 1.0/(double)sumN;
  if(sumN >= cols * (cols - 1))
    for(i=0; i<signature_len; i++) {
      signature[i+len]=(double)N[i] * wy;
      if(N[i]==0)
        signature[i]=0;
      else {
        wx = 1.0/(double)N[i];
        signature[i] = sqrt(wx*(ssx[i] - wx*sx[i]*sx[i] + ssy[i] - wx*sy[i]*sy[i]));
      }
    }

  free(sx);
  free(ssx);
  free(sy);
  free(ssy);
  free(N);

  if(sumN < cols * (cols - 1))
    return 0;
  else {
    return 1;
  }
}

int jcov_len(EZGDAL_LAYER **layers, int num_of_layers) {
    if(layers[0]->stats==NULL) return -1;
    if(layers[0]->stats->map_max_val==0) return -1;

    return ((layers[0]->stats->map_max_val+1)*(layers[0]->stats->map_max_val+2));
}

