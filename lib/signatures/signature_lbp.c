#include <math.h>
#include <ezgdal.h>
#include <assert.h>
#include <stdarg.h>


int lbp_recode[256] = {0,1,1,2,1,3,2,4,1,5,
                       3,6,2,7,4,8,1,9,5,10,
                       3,11,6,12,2,13,7,14,4,15,
                       8,16,1,5,9,13,5,17,10,18,
                       3,17,11,19,6,20,12,21,2,10,
                       13,22,7,23,14,24,4,18,15,25,
                       8,26,16,27,1,3,5,7,9,11,
                       13,15,5,17,17,20,10,23,18,26,
                       3,11,17,23,11,28,19,29,6,19,
                       20,30,12,29,21,31,2,6,10,14,
                       13,19,22,25,7,20,23,30,14,30,
                       24,32,4,12,18,24,15,29,25,33,
                       8,21,26,32,16,31,27,34,1,2,
                       3,4,5,6,7,8,9,10,11,12,
                       13,14,15,16,5,13,17,18,17,19,
                       20,21,10,22,23,24,18,25,26,27,
                       3,7,11,15,17,20,23,26,11,23,
                       28,29,19,30,29,31,6,14,19,25,
                       20,30,30,32,12,24,29,33,21,32,
                       31,34,2,4,6,8,10,12,14,16,
                       13,18,19,21,22,24,25,27,7,15,
                       20,26,23,29,30,31,14,25,30,32,
                       24,33,32,34,4,8,12,16,18,21,
                       24,27,15,26,29,31,25,32,33,34,
                       8,16,21,27,26,31,32,34,16,27,
                       31,34,27,34,34,35};

int direction[8][2] = {{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{1,0}};

int lbp(EZGDAL_LAYER *l, EZGDAL_FRAME *f, int r, int c) {
  double v,v1;
  int i, b;

  v = f->buffer[r][c];
  if(ezgdal_is_null(l,v)) return -1;

  b=0;
  for(i=0; i<8; i++) {
    v1 = f->buffer[r+direction[i][0]][c+direction[i][1]];
    if(ezgdal_is_null(l,v1)) return -1;
    if(v>v1) b=b | 1<<i;
  }

  return lbp_recode[b];
}

int local_binary_pattern(EZGDAL_FRAME **frames, int num_of_frames, double *signature, int signature_len, ...) {
  int i, r, c, cat, rows, cols, N;
  EZGDAL_LAYER *l;
  EZGDAL_FRAME *f;

  for(i=0; i<signature_len; i++)
    signature[i]=0.0;


  f = frames[0];
  N = 0;
  l = f->owner.stripe->layer;
  cols = f->rows;
  rows = f->cols;

  for(r=1; r<rows-1; r++)
    for(c=1; c<cols-1; c++) {

      cat = lbp(l,f,r,c);
      if(cat!=-1) {
        signature[cat]+=1.0;
        N++;
      }
    }

  if(2 * N < (rows-2) * (cols-2))
    return 0;
  else
    return 1;

}

int local_binary_pattern_len(EZGDAL_LAYER **layers, int num_of_layers) {
    return 36;
}

