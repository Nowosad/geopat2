/****************************************************************************
 *
 * MODULE:	Entropy of first layer
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


#define SIGNATURE_H_MAX_N 409600


/**
 * The signature contains 3 columns:
 * 1 - entropy
 * 2 - no of elementst in segment
 * 3 - size of segment
 */



typedef struct  {
  int id;
  int count;
} SIGNATURE_H_ELEMENT;

void insert_cell(SIGNATURE_H_ELEMENT *elements, int *N, int id) {
  int i;
  
  for(i=0; i<(*N); i++)
    if(elements[i].id == id) {
      elements[i].count++;
      return;
    }

  assert(*N<SIGNATURE_H_MAX_N);

  elements[*N].id = id;
  elements[*N].count = 1;
  (*N)++;
}


int H(EZGDAL_FRAME **frames, int num_of_frames, double *signature, int signature_len, ...) {
  int N, i, N_elements;
  int r, c, rows, cols;
  double **buf, sum, x, w;
  SIGNATURE_H_ELEMENT *elements = calloc(SIGNATURE_H_MAX_N,sizeof(SIGNATURE_H_ELEMENT));
  EZGDAL_LAYER *l = frames[0]->owner.frameset->layer;
  
  
  N = 0;
  N_elements = 0;
  rows = frames[0]->rows;
  cols = frames[0]->cols;
  buf = frames[0]->buffer;
  
  for(r=0; r<rows; r++)
    for(c=0; c<cols; c++) 
      if(!ezgdal_is_null(l,buf[r][c])) {
        insert_cell(elements, &N_elements, (int)buf[r][c]);
        N++;
      }

  signature[1] = N_elements;
  signature[2] = N; 

/*
  if(2*N < cols * rows) {
    free(elements);
  signature[3] = 1;
    return 0;
  }
*/
  w = (double)1.0/(double)N;
  sum = 0.0;
  for(i=0; i<N_elements; i++) {
    x = w*elements[i].count;
    sum -= x*log(x);
  }
  sum /= log(2);
  
  signature[0] = sum;
  
  free(elements);
  return 1;
}

int H_len(EZGDAL_LAYER **layers, int num_of_layers) {
    return 3;
}