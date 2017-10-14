/****************************************************************************
 *
 * MODULE:	Cartesian product
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

int cartesianproduct(EZGDAL_FRAME **frames, int num_of_frames, double *signature, int signature_len, ...) {
  int i, r, c, mx, ct, cat, N;
  EZGDAL_LAYER *l;
  EZGDAL_FRAME *f;
  double v;
  int cell_null = FALSE;

  for(i=0; i<signature_len; i++)
    signature[i]=0.0;


  f = frames[0];
  N = 0;
  l = f->owner.stripe->layer;

  for(r=0; r<f->rows; r++) 
    for(c=0; c<f->cols; c++) {

      v = f->buffer[r][c];
      cell_null = ezgdal_is_null(l,v);

      if(!cell_null) {

        cat = ezgdal_get_value_index(l, v);
        assert(cat>=0);
        for(i=num_of_frames-1; i>0; i--) {

          v = frames[i]->buffer[r][c];
          cell_null = cell_null || ezgdal_is_null(l,v);
          if(!cell_null) {
            l = frames[i]->owner.stripe->layer;
            mx = l->stats->map_max_val;
            cat*=mx;
            ct = ezgdal_get_value_index(l, v);
            assert(ct>=0);
            cat+=ct;
          }
        }

        assert(cat<signature_len);
        signature[cat]+=1.0;
        N++;
      }
    }

  if(2 * N < f->rows * f->cols)
    return 0;
  else
    return 1;
}

int cartesianproduct_len(EZGDAL_LAYER **layers, int num_of_layers) {
  int i, l=1;

  for(i=0; i<num_of_layers; i++) {
    if(layers[i]->stats==NULL) return -1;
    if(layers[i]->stats->map_max_val==0) return -1;
    l*=layers[i]->stats->map_max_val+1;
  }

  return l;
}

