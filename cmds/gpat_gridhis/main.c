/****************************************************************************
 *
 * PROGRAM:	gpat_gridhis - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel, Jakub Nowosad
 * PURPOSE:	program for creating a grids of motifels;
 *		functionality based on p.sig.grid from
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
#include <omp.h>

#include "../../lib/ezGDAL/ezgdal.h"
#include "../../lib/SML/sml.h"

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

    int ninputs,i,r,c;

    int size_val = 150;
    int shift_val = 100;
    int level_val;
    signature_func *sign_func = get_signature("cooc");
    signature_len_func *sign_len_func = get_signature_len("cooc");
    normalization_func *norm_func = get_normalization_method("pdf");

    struct arg_str  *inp   = arg_str1("i","input","<file_name>","name of input file (GeoTIFF)");
    struct arg_str  *out   = arg_str1("o","output","<file_name>","name of output file (GRID)");
    struct arg_str  *sig   = arg_str0("s","signature","<signature_name>","motifel's signature (use -l to list all signatures, default: 'cooc')");
    struct arg_int  *lvl   = arg_int0(NULL,"level","<n>","full decomposition level (default: 0, auto)");
    struct arg_int  *size  = arg_int0("z","size","<n>","motifel size in cells (default: 150)");
    struct arg_int  *shift = arg_int0("f","shift","<n>","shift of motifels (default: 100)");
    struct arg_str  *norm  = arg_str0("n","normalization","<normalization_name>","signature normalization method (use -l to list all methods, default: 'pdf')");
    struct arg_lit  *list  = arg_lit0("l",NULL,"list all signatures and normalization methods");
    struct arg_int  *th    = arg_int0("t",NULL,"<n>","number of threads (default: 1)");
    struct arg_lit  *help  = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end   = arg_end(20);
    void* argtable[] = {inp,out,sig,lvl,size,shift,norm,list,th,help,end};

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
        printf("List of available signatures:\n");
        char *list = list_all_signatures();
        printf("\n%s\n",list);
        free(list);
        exit(0);
      }
    }
    if(strcmp(sig->sval[0],"prod")!=0) {
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
        printf("List of available local normalization methods:\n");
        char *list = list_all_normalization_methods();
        printf("\n%s\n",list);
        free(list);
        exit(0);
      }
    }

    if(size->count > 0) {
      size_val = size->ival[0];
      if(size_val<10) {
        printf("\nSize can not be less then 10\n\n");
        exit(0);
      }
    }

    if(sig->count>0 && strcmp(sig->sval[0],"fdec")==0) {
      if((size_val != 0) && ((size_val & (~size_val + 1)) != size_val)) {
        printf("\nFor the full decomposition, size has to be a power of two.\n\n");
        exit(0);
      }
      
      // calculate the full decomposition level
      int max_level_val = log2(size_val);
        
      if(lvl->count>0) {
        level_val = lvl->ival[0];
        if(level_val<0 || level_val>max_level_val) {
          printf("\nFor the full decomposition, 2^level cannot be greater than the size.\nThe 'level' parameter is corrected by program.\n\n");
          level_val = max_level_val;
        } 
      } else {
          // if level is not set
          level_val = max_level_val;
      }
    }

    if(shift->count > 0) {
      shift_val = shift->ival[0];
      if(shift_val<5) {
        printf("\n'shift' can not be less then 5\n\n");
        exit(0);
      }
    }

    if(shift_val>size_val) {
      printf("\n'shift' can not be greater then 'size' parameter\n\n");
      exit(0);
    }

    /* set number of threads */
    if (th->count > 0) 
      omp_set_num_threads(th->ival[0]);
    else
      omp_set_num_threads(1);


    for(i=0; i<inp->count; i++)
      if(!ezgdal_file_exists((char *)(inp->sval[i]))) {
        printf("\nFile '%s' does not exists!\n\n", inp->sval[i]);
        usage(argv[0],argtable);
      }



    EZGDAL_LAYER **input_layers;

    input_layers = (EZGDAL_LAYER **)malloc(ninputs*sizeof(EZGDAL_LAYER *));

    for(i=0; i<ninputs; i++) {
      input_layers[i] = ezgdal_open_layer((char *)(inp->sval[i]));
      if(input_layers[i]==NULL) {
        printf("\nCannot open file: '%s'\n\n", inp->sval[i]);
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

    for(i=0; i<ninputs; i++) {
      ezgdal_create_stripe(input_layers[i],0,size_val);
      ezgdal_create_all_frames(input_layers[i]->stripe,0,shift_val);
    }

    printf("Calculating statistics...     "); fflush(stdout);
#pragma omp parallel for
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
    dims[0] = sign_len_func(input_layers, ninputs, level_val);
    if(dims[0]<0) {
      printf("\nSignature length cannot be determined!\n\n");
      usage(argv[0],argtable);
    }
    printf("Signature length: %d\n",dims[0]); fflush(stdout);

    double *at = ezgdal_layer_get_at(input_layers[0]);
    char *wkt = ezgdal_layer_get_wkt(input_layers[0]);

    SML_CELL_TYPE *cell_type = sml_create_cell_type(SML_DOUBLE,1,dims);
    SML_WINDOW *window = sml_create_window(	
                                    1+(input_layers[0]->rows-size_val)/shift_val,
                                    1+(input_layers[0]->cols-size_val)/shift_val,
                                    at[0]+at[1]*floor((size_val-shift_val)/2),
                                    at[1]*shift_val,
                                    at[2],
                                    at[3]+at[5]*floor((size_val-shift_val)/2),
                                    at[4],
                                    at[5]*shift_val,
                                    wkt
                                  );

    SML_DATA_HEADER *dh = sml_create_layer((char *)(out->sval[0]), cell_type, window);
    sml_set_layer_description(dh, argv, argc);

    void *buf = sml_create_cell_row_buffer(dh);
    printf("Calculating grid of signatures...     "); fflush(stdout);
    for(r=0; r<dh->file_win->rows; r++) {
//printf("r: %d/%d\n",r,dh->file_win->rows); fflush(stdout);
      ezgdal_show_progress(stdout,r,dh->file_win->rows);
//printf("ninputs: %d\n",ninputs); fflush(stdout);
      for(i=0; i<ninputs; i++)
        ezgdal_load_stripe_data(input_layers[i]->stripe,r*shift_val);
//printf("data read\n"); fflush(stdout);

//#pragma omp parallel for 
      for(c=0; c<dh->file_win->cols; c++) {
        int i;

        EZGDAL_FRAME **frames = malloc(ninputs*sizeof(EZGDAL_FRAME *));
        for(i=0; i<ninputs; i++)
          frames[i] = &(input_layers[i]->stripe->frame[c]);

        void *cell=sml_get_cell_pointer(dh,buf,c);
        double *cell_data = sml_get_cell_data(cell);
        if(sign_func(frames, ninputs, cell_data, dims[0], level_val)==1) {

          sml_set_cell_not_null(cell);
          if(norm_func(cell_data, dims[0])!=0) {
            printf("\nNormalization error!\n\n");
            sml_set_cell_null(cell);
          }

        } else
          sml_set_cell_null(cell);

        free(frames);
      }
      sml_write_next_row_to_layer(dh,buf);
    }
    ezgdal_show_progress(stdout,100,100);

    sml_close_layer(dh);
    free(buf);

    for(i=0; i<ninputs; i++) 
      ezgdal_close_layer(input_layers[i]);
    free(input_layers);

return 0;

}

