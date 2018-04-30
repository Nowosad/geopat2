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
#include "../../lib/ezGDAL/ezgdal.h"
#include <assert.h>
#include <stdarg.h>

int triangular_index(int r, int c) {
  r++; c++;
  if(c<=r)
    return (r-1)*r/2+c-1;
  else
    return (c-1)*c/2+r-1;

//  return (n*c1)+c2-((c1*(c1+1))/2);
}


int coocurrence(EZGDAL_FRAME **frames, int num_of_frames, double *signature, int signature_len, ...) {
  int i, r, c, cat1, cat2, N, cols, rows;
  EZGDAL_LAYER *l;
  EZGDAL_FRAME *f;
  double v;

  for(i=0; i<signature_len; i++)
    signature[i]=0.0;


  f = frames[0];
  N = 0;
  l = f->owner.stripe->layer;
  cols = f->cols;
  rows = f->rows;

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

          i = triangular_index(cat1,cat2);
          signature[i]+=1.0;
          N++;
        }

        v = f->buffer[r][c+1];
        if(!ezgdal_is_null(l,v)) {

          cat2 = ezgdal_get_value_index(l, v);
          assert(cat2>=0);

          i = triangular_index(cat1,cat2);
          assert(i>=0 && i<signature_len);
          signature[i]+=1.0;
          N++;
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

        i = triangular_index(cat1,cat2);
        signature[i]+=1.0;
        N++;
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

        i = triangular_index(cat1,cat2);
        signature[i]+=1.0;
        N++;
      }
    }
  }

  if(N < cols * (cols - 1))
    return 0;
  else
    return 1;

}

int coocurrence_len(EZGDAL_LAYER **layers, int num_of_layers) {
    if(layers[0]->stats==NULL) return -1;
    if(layers[0]->stats->map_max_val==0) return -1;

    return ((layers[0]->stats->map_max_val+1)*(layers[0]->stats->map_max_val+2))/2;
}

