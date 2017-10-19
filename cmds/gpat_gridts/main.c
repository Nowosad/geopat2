/****************************************************************************
 *
 * PROGRAM:	gpat_gridts - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	program for creating a grids of time series;
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
#include <math.h>
#include <omp.h>

#include <sml.h>
#include <ezgdal.h>


#include "../../lib/argtable/argtable3.h"
#include "../../lib/tools/libtools.h"


void usage(char *progname, void *argtable) {
      printf("\nUsage:\n\t%s", progname);
      arg_print_syntax(stdout,argtable,"\n");
      printf("\n");
      arg_print_glossary_gnu(stdout,argtable);
      printf("\n");
      exit(0);
}


int main(int argc, char **argv) {

    int i, row, col;

    int dim_val = 1;

    struct arg_str  *inp  = arg_strn("i","input","<file_name>",1,9999,"name of input file(s) (GeoTIFF)");
    struct arg_str  *out  = arg_str1("o","output","<file_name>","name of output file (GRID)");
    struct arg_int  *dim  = arg_int0("d","dimension","<n>","dimension of vector that describes time series element (default: 1)");
    struct arg_lit  *norm = arg_lit0("n","normalize","normalize each vector coordinate to [0.0, 1.0] (default: no)");
    struct arg_lit  *help = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end  = arg_end(20);
    void* argtable[] = {inp,out,dim,norm,help,end};

    int nerrors = arg_parse(argc,argv,argtable);

    if (help->count > 0) 
      usage(argv[0],argtable);

    if (nerrors > 0) {
      arg_print_errors(stdout,end,argv[0]);
      usage(argv[0],argtable);
    }

    if(dim->count>0) {
      dim_val = dim->ival[0];
      if(inp->count % dim_val != 0) {
        printf("\nNumber of input files does not fit to vector dimension!\n\n");
        usage(argv[0],argtable);
      }
    }

    for(i=0; i<inp->count; i++)
      if(!ezgdal_file_exists((char *)(inp->sval[i]))) {
        printf("\nFile '%s' does not exists!\n\n", inp->sval[i]);
        usage(argv[0],argtable);
      }

    EZGDAL_LAYER **input_layers;

    input_layers = (EZGDAL_LAYER **)malloc(inp->count*sizeof(EZGDAL_LAYER *));

    for(i=0; i<inp->count; i++) {
      input_layers[i] = ezgdal_open_layer((char *)(inp->sval[i]));
      if(input_layers[i]==NULL) {
        printf("\nCannot open file: '%s'\n\n", inp->sval[i]);
        usage(argv[0],argtable);
      }
    }

    if(!ezgdal_is_projection_ok(input_layers, inp->count)) {
      printf("\nInput files have various projections!\n\n");
      usage(argv[0],argtable);
    }

    if(!ezgdal_is_bbox_ok(input_layers,inp->count)) {
      printf("\nInput files have various extents and/or resolutions!\n\n");
      usage(argv[0],argtable);
    }

    if(norm->count>0) 
      for(i=0; i<inp->count; i++) {
        ezgdal_calc_layer_stats(input_layers[i]);
        if(norm->count>0 && input_layers[i]->stats->min==input_layers[i]->stats->max) {
          printf("\nInput layer (%d) contains only one value!\n\n", i);
          usage(argv[0],argtable);
        }
      }


    int *dims = (int *)malloc(2*sizeof(int));
    dims[0] = dim_val;
    dims[1] = inp->count / dim_val;

    SML_CELL_TYPE *cell_type = sml_create_cell_type(SML_DOUBLE,2,dims);
    double *at = ezgdal_layer_get_at(input_layers[0]);
    char *wkt = ezgdal_layer_get_wkt(input_layers[0]);
    SML_WINDOW *window = sml_create_window(
                                    input_layers[0]->rows,
                                    input_layers[0]->cols,
                                    at[0],at[1],at[2],at[3],at[4],at[5],
                                    wkt
                                  );
    SML_DATA_HEADER *dh = sml_create_layer((char *)(out->sval[0]), cell_type, window);

    if(dh!=NULL) {
      void *buf = sml_create_cell_row_buffer(dh);

      /* import data */
      ezgdal_show_progress(stdout,0,dh->file_win->rows);
      for(row=0; row<dh->file_win->rows; row++) {

        ezgdal_show_progress(stdout,row,dh->file_win->rows);
        for(i=0; i<inp->count; i++) 
            ezgdal_read_buffer(input_layers[i], row);

        for(col=0; col<dh->file_win->cols; col++) {
          void *cell = sml_get_cell_pointer(dh, buf, col);
          sml_set_cell_not_null(cell);

          for(i=0; i<inp->count; i++) {

            if(sml_is_cell_null(cell)) 
               continue;

            double v;
            v = input_layers[i]->buffer[col];
            if(ezgdal_is_null(input_layers[i],v))
              sml_set_cell_null(cell);
            else {
              if(norm->count>0)
                v = (v - input_layers[i]->stats->min)/(input_layers[i]->stats->max - input_layers[i]->stats->min);

              sml_set_cell_val_dbl(dh, v, cell, i);
            }
          }
        }
        sml_write_next_row_to_layer(dh,buf);
      }
      ezgdal_show_progress(stdout,100,100);
      sml_set_layer_description(dh,argv,argc);

    free(buf);

}

    sml_close_layer(dh);

    for(i=0; i<inp->count; i++) 
      ezgdal_close_layer(input_layers[i]);
    free(input_layers);


return 0;

}

