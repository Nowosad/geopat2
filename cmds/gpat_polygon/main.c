/****************************************************************************
 *
 * PROGRAM:	gpat_polygon - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	program for calculating signatures of polygons;
 *		functionality based on p.sig.polygon from
 *		GRASS GeoPAT by Jasiewicz, Netzel, Stepinski
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
#include <float.h>
#include <omp.h>

#include <sml.h>
#include <ezgdal.h>


#include "../../lib/argtable/argtable3.h"
#include "../../lib/signatures/signatures.h"
#include "../../lib/normalization/methods.h"
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

    int ninputs,i,j,r,c;

    signature_func *sign_func = get_signature("cooc");
    signature_len_func *sign_len_func = get_signature_len("cooc");
    normalization_func *norm_func = get_normalization_method("pdf");

    struct arg_str  *inp   = arg_str1("i","input","<file_name>","name of input file (GeoTIFF)");
    struct arg_str  *seg   = arg_str1("e","segments","<file_name>","name of input file (GeoTIFF, int)");
    struct arg_str  *out   = arg_str1("o","output","<file_name>","name of output file (TXT)");
    struct arg_str  *sig   = arg_str0("s","signature","<signature_name>","signature method (use -l to list all methods, default: 'cooc')");
    struct arg_str  *norm  = arg_str0("n","normalization","<normalization_name>","signature normalization method (use -l to list all methods, default: 'pdf')");
    struct arg_lit  *list  = arg_lit0("l",NULL,"list all signatures and normalization methods");
    struct arg_int  *th    = arg_int0("t",NULL,"<n>","number of threads (default: 1)");
    struct arg_lit  *help  = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end   = arg_end(20);
    void* argtable[] = {inp,seg,out,sig,norm,list,th,help,end};

    int nerrors = arg_parse(argc,argv,argtable);

    if (help->count > 0) 
      usage(argv[0],argtable);

    /* list all measures */
    if(list->count > 0) {
      char *list = list_all_signatures();
      printf("\nList of signatures:\n\n%s\n",list);
      free(list);
      list = list_all_normalization_methods();
      printf("\nList of local normalization methods:\n\n%s\n",list);
      free(list);
      exit(0);
    }
    
    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0) {
      /* Display the error details contained in the arg_end struct.*/
      arg_print_errors(stdout,end,argv[0]);
      usage(argv[0],argtable);
    }

    if(sig->count > 0) {
      sign_func = get_signature((char *)(sig->sval[0]));
      sign_len_func = get_signature_len((char *)(sig->sval[0]));
      /* signature not found */
      if(sign_func==NULL || sign_len_func==NULL) {
        printf("\nWrong signature name: %s\n\n",sig->sval[0]);
        printf("list of avalable signatures:\n");
        char *list = list_all_signatures();
        printf("\n%s\n",list);
        free(list);
        exit(0);
      }
    }
    if(strcmp(sig->sval[0],"cpr")!=0) {
      if(inp->count>1)
        printf("\nSignatures will be calculated for [%s] only!\n\n",inp->sval[0]);
      ninputs=1;
    } else {
      ninputs=inp->count;
    }

    if(norm->count > 0) {
      norm_func = get_normalization_method((char *)(norm->sval[0]));
      /* signature not found */
      if(norm_func==NULL) {
        printf("\nWrong signature name: %s\n\n",norm->sval[0]);
        printf("list of avalable local normalization methods:\n");
        char *list = list_all_normalization_methods();
        printf("\n%s\n",list);
        free(list);
        exit(0);
      }
    }

    /* set number of threads */
    if (th->count > 0) 
      omp_set_num_threads(th->ival[0]);
    else
      omp_set_num_threads(1);


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

    input_layers = (EZGDAL_LAYER **)malloc(ninputs*sizeof(EZGDAL_LAYER *));

    input_layers[0] = ezgdal_open_layer((char *)(seg->sval[0]));
    if(input_layers[0]==NULL) {
      printf("\nCan not open file: '%s'\n\n", seg->sval[0]);
      usage(argv[0],argtable);
    }

    for(i=1; i<ninputs; i++) {
      input_layers[i] = ezgdal_open_layer((char *)(inp->sval[i-1]));
      if(input_layers[i]==NULL) {
        printf("\nCan not open file: '%s'\n\n", inp->sval[i-1]);
        usage(argv[0],argtable);
      }
    }


    if(!ezgdal_is_projection_ok(input_layers,ninputs)) {
      printf("\nInput files have various projections!\n\n");
      usage(argv[0],argtable);
    }

    if(!ezgdal_is_bbox_ok(input_layers,ninputs)) {
      printf("\nInput files have various extent and/or resolution!\n\n");
      usage(argv[0],argtable);
    }

    printf("Calculating statistics ... "); fflush(stdout);
    for(i=0; i<ninputs; i++) {
      ezgdal_calc_layer_stats(input_layers[i]);
      double min = input_layers[i]->stats->min;
      min = min - 0.5;
      double max = input_layers[i]->stats->max;
      max = max + 0.5;
      int N = (int)floor(fabs(max-min));
      ezgdal_calc_value_map(input_layers[i],min,max,N);
    }
    printf("OK\n"); fflush(stdout);


    int *dims = (int *)malloc(sizeof(int));
    dims[0] = sign_len_func(input_layers+1, ninputs-1);
    if(dims[0]<0) {
      printf("\nSignature lenght can not be determined!\n\n");
      usage(argv[0],argtable);
    }
    printf("Signature length: %d\n",dims[0]); fflush(stdout);

    FILE *file = fopen(out->sval[0],"w");


//////////////////////////////////////////////////////////////////////////////

    printf("Calculating extents of polygons ...     "); fflush(stdout);

    EZGDAL_LAYER *lCats = input_layers[0];
    int nCats = lCats->stats->map_max_val+1;
    EZGDAL_FRAMESET *frmsetCats = ezgdal_create_frameset_with_size(lCats,nCats);
    for(i=0; i<nCats; i++) {
      EZGDAL_FRAME *frm = ezgdal_add_frameset_frame(frmsetCats,lCats->cols,0,lCats->rows,0);
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

    printf("Calculating signatures of polygons ...     "); fflush(stdout);

    for(i=1; i<ninputs; i++) {
      EZGDAL_FRAMESET *fs = ezgdal_create_frameset_with_size(input_layers[i],nCats);
      for(j=0; j<nCats; j++)
        ezgdal_add_frameset_frame(fs,
                                frmsetCats->frame[j]->col1,
                                frmsetCats->frame[j]->col2,
                                frmsetCats->frame[j]->row1,
                                frmsetCats->frame[j]->row2
                               );
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

      sign_func(frames ,ninputs-1, result, *dims);
      norm_func(result,*dims);

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

