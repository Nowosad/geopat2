/****************************************************************************
 *
 * PROGRAM:	gpat_search - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	program for calculating similarity layer;
 *		functionality based on p.sim.search from
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


void calc_simil_layer(SML_DATA_HEADER *dh, char *fname, double *refbuf, char *dtype, PALETTE *pal, int *nodata, distance_func *func) {
  int r,c,rows,cols,size;
  void *rowbuf; //, *cell;

  EZGDAL_LAYER *l;

  double a = 1.0;
  double b = 0.0;

  if(pal!=NULL) {
    a = pal->A;
    b = pal->B;
  }

  rows = dh->file_win->rows;
  cols = dh->file_win->cols;
  size = dh->cell_N_elements;
  
  rewind(dh->f);

  l = ezgdal_create_layer(fname,
                        dh->file_win->proj,
                        dtype,
                        dh->file_win->at,
                        rows,
                        cols,
                        nodata);

  rowbuf = sml_create_cell_row_buffer(dh);
 
  ezgdal_show_progress(stdout,0,rows);
  for(r=0; r<rows; r++) {
    ezgdal_show_progress(stdout,r,rows);
    sml_read_row_from_layer(dh,rowbuf,r);

#pragma omp parallel for
    for(c=0; c<cols; c++) {
      void *cell = sml_get_cell_pointer(dh,rowbuf,c);
      if(sml_is_cell_null(cell)) {
        ezgdal_set_null(l,&(l->buffer[c]));
      } else {
        double *buf[2];
        buf[0] = refbuf;
        buf[1] = sml_get_cell_data(cell);
        l->buffer[c] = a*(1.0-func(buf,2,size,1,&size))+b;
      }
    }
    ezgdal_write_buffer(l,r);
  }
  ezgdal_show_progress(stdout,100,100);

  if(strcmp(dtype,"Byte")==0)
    ezgdal_set_palette255(l,pal->pal,pal->ncolors);

  ezgdal_close_layer(l);
  
  free(rowbuf);

}


int main(int argc, char **argv) {

    char *list_dist;
    SML_DATA_HEADER *dh;
    SML_CELL_TYPE *ct;

    int size, i;
    double *refbuf;
    double x,y;
    char desc[MAX_DESC_LEN];
    char *fname;
    FILE *f;

    int *nodata = NULL;
    int _nodata = -9999;
    char *dtype = "Float64";
    distance_func *func = get_distance("jsd");
    PALETTE *palette = NULL;

    struct arg_str  *inp   = arg_str1("i","input","<file_name>","name of input file (GRID)");
    struct arg_str  *out   = arg_str0("o","output","<file_name>","name of output file (TIFF)");
    struct arg_str  *ref   = arg_str1("r","reference","<file_name>","reference data to calculate similarity (TXT)");
    struct arg_str  *mes   = arg_str0("m","measure","<measure_name>","similarity measure (use -l to list all measures; default 'jsd')");
    struct arg_lit  *mesl  = arg_lit0("l","list_measures","list all measures");
    struct arg_str  *pal   = arg_str0("p","palette","<file_name>","name of the file with colors definition (CSV)");
    struct arg_str  *type  = arg_str0(NULL,"type","Byte/....","output data type (default: Float64)");
    struct arg_int  *nodat = arg_int0("n","no_data","<n>","output NO DATA value (default: none)");
    struct arg_int  *th    = arg_int0("t",NULL,"<n>","number of threads (default: 1)");
    struct arg_lit  *help  = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end   = arg_end(20);
    void* argtable[] = {inp,out,ref,mes,mesl,type,pal,nodat,th,help,end};

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

    if(inp->count>0 && !ezgdal_file_exists((char *)(inp->sval[0]))) {
      printf("\nFile '%s' does not exists!\n\n", inp->sval[0]);
      usage(argv[0],argtable);
    }

    if(pal->count>0 && !ezgdal_file_exists((char *)(pal->sval[0]))) {
      printf("\nFile '%s' does not exists!\n\n", pal->sval[0]);
      usage(argv[0],argtable);
    }

    if(ref->count>0 && !ezgdal_file_exists((char *)(ref->sval[0]))) {
      printf("\nFile '%s' does not exists!\n\n", ref->sval[0]);
      usage(argv[0],argtable);
    }

    if(type->count>0) {
      dtype = (char *)(type->sval[0]);
      if(strcmp(dtype,"Byte")==0) {
        if(pal->count>0)
          palette = read_palette(pal->sval[0]);
        else 
          palette = create_default_palette();
      }
    }

    nodata=&_nodata;

    if(nodat->count>0) 
      _nodata = nodat->ival[0];

      dh = sml_open_layer((char *)(inp->sval[0]));
      ct = (SML_CELL_TYPE *)calloc(1,sizeof(SML_CELL_TYPE));

      size = dh->cell_N_elements;
      refbuf = malloc(size*sizeof(double));

      f = fopen(ref->sval[0],"r");
      i=1;
      while(!feof(f)) {
        if(0<=sml_read_dblbuf_txt(f,&x,&y,desc,refbuf,size,ct)) {
          fname = NULL;
          if(out->count==0)
            fname = create_fname(desc);
          else
            fname = create_fname((char *)(out->sval[0]));
          if(fname!=NULL) {
            calc_simil_layer(dh, fname, refbuf, dtype, palette, nodata, func);
            free(fname);
          }
        }
        i++;
      }

      free(refbuf);
      fclose(f);

    sml_free_cell_type(ct);
    sml_close_layer(dh);

return 0;

}

