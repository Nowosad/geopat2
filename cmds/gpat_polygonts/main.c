/****************************************************************************
 *
 * PROGRAM:	gpat_polygonts - part of GeoPAT 2
 * AUTHOR(S):	Jakub Nowosad
 * PURPOSE:	program for calculating signatures of polygons of time series;
 * COPYRIGHT:	(C) Jakub Nowosad, Space Informatics Lab,
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

#include "../../lib/ezGDAL/ezgdal.h"
#include "../../lib/SML/sml.h"

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
    struct arg_str  *seg  = arg_str1("e","segments","<file_name>","name of input file (GeoTIFF, int)");
    struct arg_str  *out  = arg_str1("o","output","<file_name>","name of output file (GRID)");
    struct arg_int  *dim  = arg_int0("d","dimension","<n>","dimension of vector that describes time series element (default: 1)");
    struct arg_lit  *norm = arg_lit0("n","normalize","normalize each vector coordinate to [0.0, 1.0] (default: no)");
    struct arg_int  *max  = arg_int0("m","max_buffer_size","<size in MB>","max size of the internal buffer for a polygon's extent, default: '4096')");
    struct arg_int  *th   = arg_int0("t",NULL,"<n>","number of threads (default: 1)");
    struct arg_lit  *help = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end  = arg_end(20);
    void* argtable[] = {inp,seg,out,dim,norm,max,th,help,end};

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
    
    if(max->count > 0) {
      if(max->ival[0]>0)
        ezgdal_frameset_max_buffer_size((unsigned long)(max->ival[0])*1048575);
    }
    
    /* set number of threads */
    if (th->count > 0) 
      omp_set_num_threads(th->ival[0]);
    else
      omp_set_num_threads(1);
    
    input_layers[0] = ezgdal_open_layer((char *)(seg->sval[0]));
    if(!ezgdal_file_exists((char *)(seg->sval[0]))) {
      printf("\nFile '%s' does not exists!\n\n", seg->sval[0]);
      usage(argv[0],argtable);
    }
    
    for(i=0; i<inp->count; i++)
      if(!ezgdal_file_exists((char *)(inp->sval[i]))) {
        printf("\nFile '%s' does not exists!\n\n", inp->sval[i]);
        usage(argv[0],argtable);
      }

    EZGDAL_LAYER **input_layers;
      
    ninputs++;
    
    input_layers = (EZGDAL_LAYER **)malloc(inp->count*sizeof(EZGDAL_LAYER *));

    for(i=1; i<=ninputs; i++) {
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
      
    /*int *dims = (int *)malloc(sizeof(int));
    dims[0] = sign_len_func(input_layers+1, ninputs-1);
    if(dims[0]<0) {
      printf("\nSignature length cannot be determined!\n\n");
      usage(argv[0],argtable);
    }
    printf("Signature length: %d\n",dims[0]); fflush(stdout);*/
    
    int *dims = (int *)malloc(2*sizeof(int));
    dims[0] = dim_val;
    dims[1] = inp->count / dim_val;
    
    FILE *file = fopen(out->sval[0],"w");
      
//////////////////////////////////////////////////////////////////////////////
      
    printf("Calculating extents of polygons...     "); fflush(stdout);
    
    EZGDAL_LAYER *lCats = input_layers[0];
    int nCats = lCats->stats->map_max_val+1;
    EZGDAL_FRAMESET *frmsetCats = ezgdal_create_frameset_with_size(lCats,nCats);
    for(i=0; i<nCats; i++) {
      EZGDAL_FRAME *frm = ezgdal_add_frameset_frame(frmsetCats,lCats->cols,0,lCats->rows,0);
      if(frm==NULL) {
        ezgdal_show_message(stderr,"Not enough RAM to proceed!");
        exit(EXIT_FAILURE);
      }
      frm->col1 = lCats->cols;
      frm->col2 = 0;
      frm->row1 = lCats->rows;
      frm->row2 = 0;
    }
    
    for(r=0; r<lCats->rows; r++) {
      ezgdal_show_progress(stdout,r,lCats->rows);
      ezgdal_read_buffer(lCats,r);
      
      for(c=0; c<lCats->cols; c++) {
        double v = lCats->buffer[c];
        if(!ezgdal_is_null(lCats,v)) {
          int idx = ezgdal_get_value_index(lCats,v);
          EZGDAL_FRAME *frm = frmsetCats->frame[idx];
          if(frm->col1 > c) frm->col1 = c;
          if(frm->col2 < c) frm->col2 = c;
          if(frm->row1 > r) frm->row1 = r;
          if(frm->row2 < r) frm->row2 = r;
        }
      }
    }
    for(i=0; i<nCats; i++) {
      EZGDAL_FRAME *frm = frmsetCats->frame[i];
      frm->cols = frm->col2 - frm->col1 + 1;
      frm->rows = frm->row2 - frm->row1 + 1;
    }
    
    ezgdal_show_progress(stdout,100,100);
    
//////////////////////////////////////////////////////////////////////////////
    
    printf("Calculating signatures of polygons...     "); fflush(stdout);
    
    for(i=1; i<ninputs; i++) {
      EZGDAL_FRAMESET *fs = ezgdal_create_frameset_with_size(input_layers[i],nCats);
      for(j=0; j<nCats; j++)
        if(ezgdal_add_frameset_frame(fs,
                                     frmsetCats->frame[j]->col1,
                                     frmsetCats->frame[j]->col2,
                                     frmsetCats->frame[j]->row1,
                                     frmsetCats->frame[j]->row2
        )==NULL) {
          ezgdal_show_message(stderr,"Not enough RAM to proceed!");
          exit(EXIT_FAILURE);
        }
        if(!input_layers[i]->is_no_data) {
          input_layers[i]->is_no_data = TRUE;
          input_layers[i]->no_data = DBL_MIN;
        }
    }
    
    double *result = (double *)malloc(dims[0]*sizeof(double));
    EZGDAL_FRAME **frames = (EZGDAL_FRAME **)malloc((ninputs-1)*sizeof(EZGDAL_FRAME *));
    
    for(i=0; i<nCats; i++) {
      ezgdal_show_progress(stdout, i, nCats);
      
      EZGDAL_FRAME *frmCats = frmsetCats->frame[i];
      ezgdal_load_frameset_frame_data(frmCats);
      
      double cat = ezgdal_get_index_value(input_layers[0],i);
      
      for(j=1; j<ninputs; j++) {
        
        EZGDAL_FRAME *frmInp = input_layers[j]->frameset->frame[i];
        frames[j-1] = frmInp;
        
        ezgdal_load_frameset_frame_data(frmInp);
        
        for(r=0; r<frmCats->rows; r++)
          for(c=0; c<frmCats->cols; c++) 
            if(ezgdal_get_value_index(input_layers[0],frmCats->buffer[r][c])!=i)
              ezgdal_set_null(input_layers[j],&(frmInp->buffer[r][c]));
            
      }
      
      double x = 0;
      double y = 0;
      int n = 0;
      
      for(r=0; r<frmCats->rows; r++)
        for(c=0; c<frmCats->cols; c++)
          if(ezgdal_get_value_index(input_layers[0],frmCats->buffer[r][c])==i) {
            x += c;
            y += r;
            n++;
          }
          c = (int)round(x/(double)n + frmCats->col1);
          r = (int)round(y/(double)n + frmCats->row1);
          x = ezgdal_cr2x(input_layers[0],c,r);
          y = ezgdal_cr2y(input_layers[0],c,r);
          char desc[64];
          sprintf(desc,"cat: %.0lf",cat);
          sml_write_dblbuf_txt(file,x,y,desc,result,*dims,1,dims);
    }
    ezgdal_show_progress(stdout, 100,100);
    
    free(result);
    free(frames);
    
    
/////////////////////////////////////////////////////////////////////////////
    
    fclose(file);
    for(i=0; i<ninputs; i++) 
      ezgdal_close_layer(input_layers[i]);
    free(input_layers);
    
    return 0;

}

