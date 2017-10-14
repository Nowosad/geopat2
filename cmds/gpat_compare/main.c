/****************************************************************************
 *
 * PROGRAM:	gpat_compare - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	program for comparing two grids of motifels;
 *		functionality based on p.sim.compare from
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

#include <sml.h>
#include <ezgdal.h>

#include "../../lib/argtable/argtable3.h"
#include "../../lib/measures/measures.h"
#include "../../lib/tools/libtools.h"

#include "palette.h"

void usage(char *progname, void *argtable) {
      printf("\nUsage:\n\t%s", progname);
      arg_print_syntax(stdout,argtable,"\n");
      printf("\n");
      arg_print_glossary_gnu(stdout,argtable);
      printf("\n");
      exit(0);
}

int is_the_same_bbox_proj(SML_DATA_HEADER *dh0, SML_DATA_HEADER *dh1) {
  return TRUE;
}


void calc_simil_layer(SML_DATA_HEADER *dh0, SML_DATA_HEADER *dh1,char *fname, char *dtype, PALETTE *pal, int *nodata, distance_func *func) {

  EZGDAL_LAYER *l;

  int r, c, rows, cols, size;
  void *rowbuf0, *rowbuf1; 

  double a = 1.0;
  double b = 0.0;

  if(pal!=NULL) {
    a = pal->A;
    b = pal->B;
  }

  rows = dh0->file_win->rows;
  cols = dh0->file_win->cols;
  size = dh0->cell_N_elements;
  
  rewind(dh0->f);
  rewind(dh1->f);

  l = ezgdal_create_layer(fname,
                        dh0->file_win->proj,
                        dtype,
                        dh0->file_win->at,
                        rows,
                        cols,
                        nodata);
                               
  rowbuf0 = sml_create_cell_row_buffer(dh0);
  rowbuf1 = sml_create_cell_row_buffer(dh1);
 
  ezgdal_show_progress(stdout,0,rows);
  for(r=0; r<rows; r++) {
    ezgdal_show_progress(stdout,r,rows);
    sml_read_row_from_layer(dh0,rowbuf0,r);
    sml_read_row_from_layer(dh1,rowbuf1,r);

#pragma omp parallel for
    for(c=0; c<cols; c++) {
      void *cell0 = sml_get_cell_pointer(dh0,rowbuf0,c);
      void *cell1 = sml_get_cell_pointer(dh1,rowbuf1,c);
      if(sml_is_cell_null(cell0) || sml_is_cell_null(cell1)) {
        ezgdal_set_null(l,&(l->buffer[c]));
      } else {
        double *buf[2];
        buf[0] = sml_get_cell_data(cell0);
        buf[1] = sml_get_cell_data(cell1);
        l->buffer[c] = a*(1.0-func(buf,2,size,dh0->cell_type->dim,dh0->cell_type->dims))+b;
      }
    }
    ezgdal_write_buffer(l,r);
  }
  ezgdal_show_progress(stdout,100,100);

  if(strcmp(dtype,"Byte")==0)
    ezgdal_set_palette255(l,pal->pal,pal->ncolors);

  free(rowbuf0);
  free(rowbuf1);

  ezgdal_close_layer(l);

}


int main(int argc, char **argv) {

    SML_DATA_HEADER *dh0, *dh1;
    char *list_dist;
	
    int *nodata = NULL;
    int _nodata = -9999;
    char *dtype = "Float64";
    distance_func *func = get_distance("jsd");
    PALETTE *palette = NULL;

    struct arg_str  *inp   = arg_strn("i","input","<file_name>",2,2,"name of input files (GRID)");
    struct arg_str  *out   = arg_str1("o","output","<file_name>","name of output file with similarity (TIFF)");
    struct arg_str  *mes   = arg_str0("m","measure","<measure_name>","similarity measure (use -l to list all measures; default 'jsd')");
    struct arg_lit  *mesl  = arg_lit0("l","list_measures","list all measures");
    struct arg_str  *pal   = arg_str0("p","palette","<file_name>","name of the file with colors definition (CSV)");
    struct arg_str  *type  = arg_str0(NULL,"type","Byte/....","output data type (default: Float64)");
    struct arg_int  *nodat = arg_int0("n","no_data","<n>","output NO DATA value (default: none)");
    struct arg_int  *th    = arg_int0("t",NULL,"<n>","number of threads (default: 1)");
    struct arg_lit  *help  = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end   = arg_end(20);
    void* argtable[] = {inp,out,type,pal,nodat,mes,mesl,th,help,end};

    int nerrors = arg_parse(argc,argv,argtable);

    if (help->count > 0) 
      usage(argv[0],argtable);

    /* list all measures */
    if(mesl->count > 0) {
      char *list = list_all_distances();
      printf("\nList of measures:\n\n%s\n",list);
      free(list);
      exit(0);
    }
    
    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0) {
      /* Display the error details contained in the arg_end struct.*/
      arg_print_errors(stdout,end,argv[0]);
      usage(argv[0],argtable);
    }

    if(mes->count > 0) {
      func = get_distance((char *)(mes->sval[0]));
      /* measure not found */
      if(func==NULL) {
        printf("\nWrong similarity measure: %s\n\n",mes->sval[0]);
        printf("list of distances:\n");
        list_dist = list_all_distances();
        printf("\n%s\n",list_dist);
        free(list_dist);
        exit(0);
      }
    }

    /* set number of threads */
    if (th->count > 0) 
      omp_set_num_threads(th->ival[0]);
    else
      omp_set_num_threads(1);

    if(inp->count!=2) {
      printf("\nTwo input files (GRID) are needed!\n\n");
      usage(argv[0],argtable);
    } else if(!ezgdal_file_exists((char *)(inp->sval[0]))) {
      printf("\nInput file '%s' does not exists!\n\n", inp->sval[0]);
      usage(argv[0],argtable);
    } else if(!ezgdal_file_exists((char *)(inp->sval[1]))) {
      printf("\nInput file '%s' does not exists!\n\n", inp->sval[1]);
      usage(argv[0],argtable);
    }


    if(type->count>0) {
      dtype = (char *)(type->sval[0]);
      if(strcmp(dtype,"Byte")==0) {
        if(pal->count>0) {
          if(!ezgdal_file_exists((char *)(pal->sval[0]))) {
            printf("\nPalette file '%s' does not exists!\n\n", pal->sval[0]);
            usage(argv[0],argtable);
          }
          palette = read_palette(pal->sval[0]);
        } else 
          palette = create_default_palette();
      } else if(pal->count>0) 
        printf("\nResult type is not 'Byte'. Palette file '%s' ignored!\n\n", pal->sval[0]);
    }

    if(nodat->count>0) {
      _nodata = nodat->ival[0];
      nodata=&_nodata;
    }
   
    dh0 = sml_open_layer((char *)(inp->sval[0]));
    dh1 = sml_open_layer((char *)(inp->sval[1]));

    if(!is_the_same_bbox_proj(dh0,dh1)) {
      printf("\nInput files have different pojection or boundig box!\n\n");
      sml_close_layer(dh0);
      sml_close_layer(dh1);
      usage(argv[0],argtable);
    }

    calc_simil_layer(dh0, dh1, (char *)(out->sval[0]), dtype, palette, nodata, func);

    sml_close_layer(dh0);
    sml_close_layer(dh1);

return 0;

}

