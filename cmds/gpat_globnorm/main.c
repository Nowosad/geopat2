/****************************************************************************
 *
 * PROGRAM:	gpat_globnorm - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	program for global normalization of vector elements;
 * COPYRIGHT:	(C) Pawel Netzel, Space Informatics Lab,
 *		University of Cincinnati
 *              http://sil.uc.edu
 *
 *		This program is free software under 
 *		the GNU General Public License (>=v3). 
 *		https://www.gnu.org/licenses/gpl-3.0.en.html
 *
 *****************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <omp.h>

#include <sml.h>
#include <ezgdal.h>

#include "../../lib/argtable/argtable3.h"
#include "../../lib/tools/libtools.h"

typedef struct {
    double *min, *max, *sum, *sum2;
    long int N;
} NORM_STATS;

NORM_STATS *create_norm_stats(SML_DATA_HEADER *dh) {
    int i,n;
    NORM_STATS *ns = (NORM_STATS *)malloc(sizeof(NORM_STATS));
    n = dh->cell_N_elements;
    ns->min = (double *)malloc(n*sizeof(double));
    ns->max = (double *)malloc(n*sizeof(double));
    ns->sum = (double *)malloc(n*sizeof(double));
    ns->sum2 = (double *)malloc(n*sizeof(double));
    ns->N = 0;
    for(i=0; i<n; i++) {
      ns->min[i] = DBL_MAX;
      ns->max[i] = DBL_MIN;
      ns->sum[i] = 0.0;
      ns->sum2[i] = 0.0;
    }
    return ns;
}

void free_norm_stats(NORM_STATS *ns) {
    free(ns->min);
    free(ns->max);
    free(ns->sum);
    free(ns->sum2);
    free(ns);
}


void usage(char *progname, void *argtable) {
      printf("\nUsage:\n\t%s", progname);
      arg_print_syntax(stdout,argtable,"\n");
      printf("\n");
      arg_print_glossary_gnu(stdout,argtable);
      printf("\n");
      exit(0);
}

void scan_stats(SML_DATA_HEADER *dh0, NORM_STATS *ns) {

  int i, r, c, rows, cols, size;
  void *rowbuf; 
  double v;

  rows = dh0->file_win->rows;
  cols = dh0->file_win->cols;
  size = dh0->cell_N_elements;
  
  rewind(dh0->f);

  rowbuf = sml_create_cell_row_buffer(dh0);
 
  ezgdal_show_progress(stdout,0,rows);
  for(r=0; r<rows; r++) {
    ezgdal_show_progress(stdout,r,rows);
    sml_read_row_from_layer(dh0,rowbuf,r);

    for(c=0; c<cols; c++) {
      void *cell = sml_get_cell_pointer(dh0,rowbuf,c);
      if(!sml_is_cell_null(cell)) {
        for(i=0; i<size; i++) {
          v = sml_get_cell_val_dbl(dh0, cell, i);
          if(ns->min[i]>v) ns->min[i] = v;
          if(ns->max[i]<v) ns->max[i] = v;
          ns->sum[i] += v;
          ns->sum2[i] += v*v;
        }
        ns->N++;
      }
    }
  }
  ezgdal_show_progress(stdout,100,100);

  free(rowbuf);
}


void normalize_01(SML_DATA_HEADER *dh0, SML_DATA_HEADER *dh1, NORM_STATS *ns) {

  int i, r, c, rows, cols, size;
  void *rowbuf0, *rowbuf1; 
  double v;

  rows = dh0->file_win->rows;
  cols = dh0->file_win->cols;
  size = dh0->cell_N_elements;
  
  rewind(dh0->f);

  rowbuf0 = sml_create_cell_row_buffer(dh0);
  rowbuf1 = sml_create_cell_row_buffer(dh1);

  for(i=0; i<size; i++)
    ns->max[i] = ns->max[i] - ns->min[i];

  ezgdal_show_progress(stdout,0,rows);
  for(r=0; r<rows; r++) {
    ezgdal_show_progress(stdout,r,rows);
    sml_read_row_from_layer(dh0,rowbuf0,r);

    for(c=0; c<cols; c++) {
      void *cell0 = sml_get_cell_pointer(dh0,rowbuf0,c);
      void *cell1 = sml_get_cell_pointer(dh1,rowbuf1,c);
      if(!sml_is_cell_null(cell0)) {
        sml_set_cell_not_null(cell1);
        for(i=0; i<size; i++) {
          v = sml_get_cell_val_dbl(dh0, cell0, i);
          if(ns->max[i]!=0.0)
            v = (v - ns->min[i])/ns->max[i];
          else 
            v = 0.0;
          sml_set_cell_val_dbl(dh1, v, cell1, i);
        }
        ns->N++;
      } else
        sml_set_cell_null(cell1);

    }

    sml_write_next_row_to_layer(dh1,rowbuf1);
  }
  ezgdal_show_progress(stdout,100,100);

  free(rowbuf0);
  free(rowbuf1);
}

void normalize_N01(SML_DATA_HEADER *dh0, SML_DATA_HEADER *dh1, NORM_STATS *ns) {

  int i, r, c, rows, cols, size;
  void *rowbuf0, *rowbuf1; 
  double v;

  rows = dh0->file_win->rows;
  cols = dh0->file_win->cols;
  size = dh0->cell_N_elements;
  
  rewind(dh0->f);

  rowbuf0 = sml_create_cell_row_buffer(dh0);
  rowbuf1 = sml_create_cell_row_buffer(dh1);

  for(i=0; i<size; i++) {
    ns->min[i] = ns->sum[i]/(double)ns->N;
    ns->max[i] = sqrt((ns->sum2[i]-(ns->sum[i] * ns->sum[i])/(double)(ns->N))/((double)(ns->N)-1.0));
  }

  ezgdal_show_progress(stdout,0,rows);
  for(r=0; r<rows; r++) {
    ezgdal_show_progress(stdout,r,rows);
    sml_read_row_from_layer(dh0,rowbuf0,r);

    for(c=0; c<cols; c++) {
      void *cell0 = sml_get_cell_pointer(dh0,rowbuf0,c);
      void *cell1 = sml_get_cell_pointer(dh1,rowbuf1,c);
      if(!sml_is_cell_null(cell0)) {
        sml_set_cell_not_null(cell1);
        for(i=0; i<size; i++) {
          v = sml_get_cell_val_dbl(dh0, cell0, i);
          if(ns->max[i]!=0.0)
            v = (v - ns->min[i])/ns->max[i];
          else 
            v = 0.0;
          sml_set_cell_val_dbl(dh1, v, cell1, i);
        }
        ns->N++;
      } else
        sml_set_cell_null(cell1);
    }

    sml_write_next_row_to_layer(dh1,rowbuf1);
  }
  ezgdal_show_progress(stdout,100,100);

  free(rowbuf0);
  free(rowbuf1);
}

#define SPLIT 57

void normalize_ind01(SML_DATA_HEADER *dh0, SML_DATA_HEADER *dh1, NORM_STATS *ns) {

  int i, r, c, rows, cols, size;
  void *rowbuf0, *rowbuf1; 
  double v, w1, w2, eps = 0.5;

  rows = dh0->file_win->rows;
  cols = dh0->file_win->cols;
  size = dh0->cell_N_elements;
  
  w1 = (1.0-eps)*size/SPLIT;
  w2 = eps*size/(size-SPLIT);
  v = sqrt(size)/sqrt(w1*w1*SPLIT+w2*w2*(size-SPLIT));
  w1 *= v;
  w2 *= v;
  w2 *= 0.01;

  rewind(dh0->f);

  rowbuf0 = sml_create_cell_row_buffer(dh0);
  rowbuf1 = sml_create_cell_row_buffer(dh1);

  for(i=0; i<SPLIT; i++)
    ns->max[i] = ns->max[i] - ns->min[i];

  ezgdal_show_progress(stdout,0,rows);
  for(r=0; r<rows; r++) {
    ezgdal_show_progress(stdout,r,rows);
    sml_read_row_from_layer(dh0,rowbuf0,r);

    for(c=0; c<cols; c++) {
      void *cell0 = sml_get_cell_pointer(dh0,rowbuf0,c);
      void *cell1 = sml_get_cell_pointer(dh1,rowbuf1,c);
      if(!sml_is_cell_null(cell0)) {
        sml_set_cell_not_null(cell1);
        for(i=0; i<SPLIT; i++) {
          v = sml_get_cell_val_dbl(dh0, cell0, i);
          if(ns->max[i]!=0.0)
            v = (v - ns->min[i])/ns->max[i];
          else 
            v = 0.0;
          v *= w1;
          sml_set_cell_val_dbl(dh1, v, cell1, i);
        }
        for(i=SPLIT; i<size; i++) {
          v = sml_get_cell_val_dbl(dh0, cell0, i);
          v *= w2;
          sml_set_cell_val_dbl(dh1, v, cell1, i);
        }
        ns->N++;
      } else
        sml_set_cell_null(cell1);
    }

    sml_write_next_row_to_layer(dh1,rowbuf1);
  }
  ezgdal_show_progress(stdout,100,100);

  free(rowbuf0);
  free(rowbuf1);
}


int main(int argc, char **argv) {

    SML_DATA_HEADER *dh0, *dh1;

    NORM_STATS *norm_stats;

    char *norm_method = "01";
    char list_methods[] = "\nList of global normalization methods:\n\n01\tnormalize coordinats to [0,1]\nN01\tnormalize coordinats to N(0,1)\nind01\tnormalize coordinats to [0,1] for 72 lanscape indices\n";

    struct arg_str  *inp   = arg_str1("i","input","<file_name>","name of input file (GRID)");
    struct arg_str  *out   = arg_str1("o","output","<file_name>","name of output file (GRID)");
    struct arg_str  *norm  = arg_str0("m","method","<method_name>","normalization method (use -l to list all methods, default: '01')");
    struct arg_lit  *norml = arg_lit0("l","list_methods","list all methods");
    struct arg_int  *th    = arg_int0("t",NULL,"<n>","number of threads (default: 1)");
    struct arg_lit  *help  = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end   = arg_end(20);
    void* argtable[] = {inp,out,norm,norml,th,help,end};

    int nerrors = arg_parse(argc,argv,argtable);

    if (help->count > 0) 
      usage(argv[0],argtable);

    /* list all methods */
    if(norml->count > 0) {
      printf(list_methods);
      exit(0);
    }

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0) {
      /* Display the error details contained in the arg_end struct.*/
      arg_print_errors(stdout,end,argv[0]);
      usage(argv[0],argtable);
    }

    if(norm->count > 0) {
      /* measure not found */
      if(strcmp((char *)(norm->sval[0]),"01")!=0 && 
         strcmp((char *)(norm->sval[0]),"N01")!=0 &&
         strcmp((char *)(norm->sval[0]),"ind01")!=0  ) {
        printf("\nWrong nomalization method: %s\n\n",norm->sval[0]);
        printf(list_methods);
        exit(0);
      } else 
        norm_method = (char *)(norm->sval[0]);
    }

    /* set number of threads */
    if (th->count > 0) 
      omp_set_num_threads(th->ival[0]);
    else
      omp_set_num_threads(1);

    if(!ezgdal_file_exists((char *)(inp->sval[0]))) {
      printf("\nInput file '%s' does not exists!\n\n", inp->sval[0]);
      usage(argv[0],argtable);
    }

    dh0 = sml_open_layer((char *)(inp->sval[0]));

    if(norm->count>0 && strcmp((char *)(norm->sval[0]),"ind01")==0 && dh0->cell_N_elements<SPLIT ) {
        printf("\nNormalization ind01 cannot be applied to a vector of size less then %d\n\n",SPLIT);
        printf(list_methods);
        exit(0);
    }

    norm_stats = create_norm_stats(dh0);
    scan_stats(dh0, norm_stats);

    dh1 = sml_create_layer((char *)(out->sval[0]),
                           sml_create_cell_type(SML_DOUBLE,dh0->cell_type->dim,dh0->cell_type->dims),
                           sml_create_window_copy(dh0->file_win));

    if(strcmp(norm_method,"01")==0)
      normalize_01(dh0, dh1, norm_stats);
    else if(strcmp(norm_method,"N01")==0)
      normalize_N01(dh0, dh1, norm_stats);
    else
      normalize_ind01(dh0, dh1, norm_stats);

    sml_set_layer_description(dh1,argv,argc);

    free_norm_stats(norm_stats);
    sml_close_layer(dh0);
    sml_close_layer(dh1);

return 0;

}

