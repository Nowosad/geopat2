/****************************************************************************
 *
 * PROGRAM:	gpat_pointshis - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	program for calculating a signature of a surrounding of a point
 *		or a set of points;
 *		functionality based on p.sig.point from
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
#include <string.h>
#include <math.h>
#include <omp.h>

#include <ezgdal.h>
#include <sml.h>

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

    char *description;
    FILE *f, *fxy;
    char desc_text[MAX_DESC_LEN];
    double coord_x,coord_y;
    int col,row, i, sign_len;
    int ninputs;

    int size_val = 150;
    int level_val = 0;
    signature_func *sign_func = get_signature("cooc");
    signature_len_func *sign_len_func = get_signature_len("cooc");
    normalization_func *norm_func = get_normalization_method("pdf");

    struct arg_str  *inp   = arg_str1("i","input","<file_name>","name of input file (GeoTIFF)");
    struct arg_str  *out   = arg_str1("o","output","<file_name>","name of output file (TXT)");
    struct arg_str  *sign  = arg_str0("s","signature","<signature_name>","motifel's signature (use -l to list all signatures, default: 'cooc')");
    struct arg_int  *lvl   = arg_int0(NULL,"level","<n>","full decomposition level (default: 0, auto)");
    struct arg_int  *size  = arg_int0("z","size","<n>","motifel size in cells (default: 150)");
    struct arg_str  *norm  = arg_str0("n","normalization","<normalization_name>","signature normalization method (use -l to list all methods, default: 'pdf')");
    struct arg_lit  *list  = arg_lit0("l",NULL,"list all signatures and normalization methods");
    struct arg_dbl  *x     = arg_dbl0("x",NULL,"<double>","x coord");
    struct arg_dbl  *y     = arg_dbl0("y",NULL,"<double>","y coord");
    struct arg_str  *desc  = arg_str0("d","description","<string>","Description of the location");
    struct arg_str  *xy    = arg_str0(NULL,"xy_file","<file_name>","name of file with coordinates (TXT)");
    struct arg_lit  *app   = arg_lit0("a","append","append results to output file");
    struct arg_lit  *help  = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end   = arg_end(20);
    void* argtable[] = {inp,out,sign,lvl,size,norm,list,x,y,desc,xy,app,help,end};


    int nerrors = arg_parse(argc,argv,argtable);

    description = "location";

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

    if(sign->count > 0) {
      sign_func = get_signature((char *)(sign->sval[0]));
      sign_len_func = get_signature_len((char *)(sign->sval[0]));
      /* signature not found */
      if(sign_func==NULL || sign_len_func==NULL) {
        printf("\nWrong signature name: %s\n\n",sign->sval[0]);
        printf("list of avalable signatures:\n");
        char *list = list_all_signatures();
        printf("\n%s\n",list);
        free(list);
        exit(0);
      }
    }
    if(strcmp(sign->sval[0],"prod")!=0) {
      if(inp->count>1)
        printf("\nSignatures will be calculated for [%s] only!\n\n",inp->sval[0]);
      ninputs=1;
    } else {
      ninputs=inp->count;
    }

    if(sign->count>0 && strcmp(sign->sval[0],"fdec")==0) {
      if((size_val != 0) && ((size_val & (~size_val + 1)) != size_val)) {
        printf("\nFor full decomposition, size has to be a power of two.\n\n");
        exit(0);
      }
    
      if(lvl->count>0) {
        level_val = lvl->ival[0];
        if(level_val<0 || (1<<(level_val-1))>size_val) {
          printf("\nFor full decomposition, 2^level cannot be greater then size.\nThe 'level' parameter is corrected by program.\n\n");
          level_val = 0;
        }
      }
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

    if(size->count > 0) {
      size_val = size->ival[0];
      if(size_val<10) {
        printf("\nSize can not be less then 10\n\n");
        exit(0);
      }
    }

    if(!((x->count>0 && y->count>0) || (xy->count>0))) {
      printf("\n%s\n\n", "User has to provide either x and y parameter or xy_file parameter.");
      usage(argv[0],argtable);
    }

    if((xy->count>0) && !ezgdal_file_exists((char *)(xy->sval[0]))) {
      printf("\nFile [%s] does not exist.\n\n", xy->sval[0]);
      usage(argv[0],argtable);
    }

    EZGDAL_FRAMESET *frameset;
    EZGDAL_FRAME *frame;
    EZGDAL_LAYER **input_layers = (EZGDAL_LAYER **)malloc(ninputs*sizeof(EZGDAL_LAYER *));
    EZGDAL_FRAME **frames = malloc(ninputs*sizeof(EZGDAL_FRAME *));

    for(i=0; i<ninputs; i++) {
      input_layers[i] = ezgdal_open_layer((char *)(inp->sval[0]));
      if(input_layers[i]==NULL) {
        printf("\nFile [%s] cannot be opened.\n\n", inp->sval[0]);
        usage(argv[0],argtable);
      }
      frameset = ezgdal_create_frameset_with_size(input_layers[i],1);
      frame = ezgdal_add_frameset_frame(frameset,0,0,0,0);
      frames[i] = frame;
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

    if(!ezgdal_is_projection_ok(input_layers,ninputs)) {
      printf("\nInput files have various projections!\n\n");
      usage(argv[0],argtable);
    }

    if(!ezgdal_is_bbox_ok(input_layers,ninputs)) {
      printf("\nInput files have various extent and/or resolution!\n\n");
      usage(argv[0],argtable);
    }

    sign_len = sign_len_func(input_layers, ninputs, level_val);
    double *sign_buf = (double *)malloc(sign_len*sizeof(double));

    if(app->count>0)
      f = fopen(out->sval[0],"a");
    else
      f = fopen(out->sval[0],"w");

    printf("Calculating signature ... "); fflush(stdout);

    if(x->count>0 && y->count>0) {
      if(desc->count>0)
        description = (char *)(desc->sval[0]);
      col = ezgdal_xy2c(input_layers[0],x->dval[0],y->dval[0]) - size_val/2;
      row = ezgdal_xy2r(input_layers[0],x->dval[0],y->dval[0]) - size_val/2;

      for(i=0; i<ninputs; i++) {
        ezgdal_frameset_set_frame(frames[i],col,col+size_val-1,row,row+size_val-1);
        ezgdal_load_frameset_frame_data(frames[i]);
      }
/*      
printf("r: %d, c: %d, n: %d, len: %d, hist: %d, mx: %d, max: %lf, min: %lf\n",row,col,ninputs,sign_len,
                                         input_layers[0]->stats->hist_N,
                                         input_layers[0]->stats->map_max_val,
                                         input_layers[0]->stats->max,
                                         input_layers[0]->stats->min
                                         ); fflush(stdout);
*/
      if(sign_func(frames, ninputs, sign_buf, sign_len, level_val)==1) {
        // calc
        if(norm_func(sign_buf, sign_len)!=0) {
          printf("\nNormalization error!\n\n");
          sign_len = 0;
        }
      } else {
        sign_len = 0;
      }
      sml_write_dblbuf_txt(f,x->dval[0],y->dval[0],description,sign_buf,sign_len,1,&sign_len);

  } else {
      fxy = fopen((char *)xy->sval[0],"r");
      int line = 1;
      while(read_xy_txt(line++,fxy,&coord_x,&coord_y,desc_text, MAX_DESC_LEN)) {


        col = ezgdal_xy2c(input_layers[0],coord_x,coord_y) - size_val/2;
        row = ezgdal_xy2r(input_layers[0],coord_x,coord_y) - size_val/2;
        for(i=0; i<ninputs; i++) {
          ezgdal_frameset_set_frame(frames[i],col,col+size_val-1,row,row+size_val-1);
          ezgdal_load_frameset_frame_data(frames[i]);
        }
      
        if(sign_func(frames, ninputs, sign_buf, sign_len, level_val)==1) {
          // calc
          if(norm_func(sign_buf, sign_len)!=0) {
            printf("\nNormalization error!\n\n");
            sign_len = 0;
          }
        } else {
          sign_len = 0;
        }
        sml_write_dblbuf_txt(f,coord_x,coord_y,desc_text,sign_buf,sign_len,1,&sign_len);

      }
      fclose(fxy);
    }

    fclose(f);
    for(i=0; i<ninputs; i++)
      ezgdal_close_layer(input_layers[i]);

printf("OK\n");

return 0;

}
