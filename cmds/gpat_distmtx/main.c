/****************************************************************************
 *
 * PROGRAM:	gpat_distmtx - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	program for calculating matrix of distances or similarities;
 *		functionality based on p.sim.distmatrix from
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

void usage(char *progname, void *argtable) {
      printf("\nUsage:\n\t%s", progname);
      arg_print_syntax(stdout,argtable,"\n");
      printf("\n");
      arg_print_glossary_gnu(stdout,argtable);
      printf("\n");
      exit(0);
}

void count_rc(FILE *f, int *rows, int *cols, SML_CELL_TYPE *ct1) {
  int c = 0;
  int size = 0;
  double *buff = (double *)malloc(65536*sizeof(double));
  char *desc = (char *)malloc(SML_DATA_DESC_LEN*sizeof(char));
  double x,y;
  SML_CELL_TYPE *ct = ct1;

  *rows=0;
  *cols=0;
  rewind(f);
  while(!feof(f)) {
    c = sml_read_dblbuf_txt(f,&x,&y,desc,buff,size,ct);
    if(c>0) {
      (*rows)++;
      if(*cols==0) {
        *cols = c;
        size = c;
        ct = (SML_CELL_TYPE *)calloc(1,sizeof(SML_CELL_TYPE));
      }
    }
  }

  free(buff);
  free(desc);
  if(*rows>0)
    sml_free_cell_type(ct);
}

typedef struct {
  double x,y;
  char desc[SML_DATA_DESC_LEN];
  double *val;
} DATA_RECORD;

typedef struct {
  int rows;
  int cols;
  int n;
  DATA_RECORD *data;
} DATA_LIST;

DATA_LIST *read_data(FILE *f, SML_CELL_TYPE *ct1) {
  int i ,n;
  DATA_RECORD *rec;
  SML_CELL_TYPE *ct;
  DATA_LIST *list = calloc(1, sizeof(DATA_LIST));;

  count_rc(f,&(list->rows),&(list->cols), ct1);

  list->data = calloc(list->rows, sizeof(DATA_RECORD));
  ct = (SML_CELL_TYPE *)calloc(1,sizeof(SML_CELL_TYPE));

  rewind(f);
  rec = list->data;
  n=0;
  while(!feof(f) && n<list->rows) {
    rec->val = malloc(list->cols * sizeof(double));
    i = sml_read_dblbuf_txt(f, &(rec->x), &(rec->y), rec->desc, rec->val, list->cols, ct);
    if(i>0) {
      n++;
      rec++;
    }
  }

  sml_free_cell_type(ct);
  return list;
}

void free_data(DATA_LIST *list) {
  int i;
  for(i=0; i<list->n; i++)
    free(list->data[i].val);
  free(list->data);
  free(list);
}

int main(int argc, char **argv) {

    DATA_LIST *list;
    SML_CELL_TYPE *ct = (SML_CELL_TYPE *)calloc(1,sizeof(SML_CELL_TYPE));
    FILE *f, *fo;
    double *buf[2];
    int i,j;
    char *list_dist;
    distance_func *func = get_distance("jsd");

    struct arg_str  *inp   = arg_str1("i","input","<file_name>","name of input file witch signatures (TXT)");
    struct arg_str  *out   = arg_str1("o","output","<file_name>","name of output file (CSV) with similarity matrix");
    struct arg_str  *mes   = arg_str0("m","measure","<measure_name>","similarity measure (use -l to list all measures; default 'jsd')");
    struct arg_lit  *mesl  = arg_lit0("l","list_measures","list all measures");
    struct arg_lit  *sim   = arg_lit0("s","similarity","output is a similarity matrix");
    struct arg_lit  *help  = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end   = arg_end(20);
    void* argtable[] = {inp,out,mes,mesl,sim,help,end};

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

    if(inp->count>0 && !ezgdal_file_exists((char *)(inp->sval[0]))) {
      printf("\nFile '%s' does not exists!\n\n", inp->sval[0]);
      usage(argv[0],argtable);
    }

    f = fopen(inp->sval[0],"r");
    list = read_data(f, ct);

    fo = fopen(out->sval[0],"w");

    ezgdal_show_progress(stdout,0,list->rows);

    if(sim->count==0)
      fprintf(fo,"\"Distance\"");
    else
      fprintf(fo,"\"Similarity\"");

    for(i=0; i<list->rows; i++) 
      fprintf(fo,",\"%s\"",list->data[i].desc);
    fprintf(fo,"\n");

    for(i=0; i<list->rows; i++) {

      ezgdal_show_progress(stdout,i,list->rows);

      buf[0] = list->data[i].val;
      buf[1] = list->data[0].val;
//printf("%s - %s\n",list->data[i].desc,list->data[0].desc);
      fprintf(fo,"\"%s\",%.15lf",
                      list->data[i].desc,
                      (sim->count==0)?func(buf,2,list->cols,ct->dim,ct->dims):1.0-func(buf,2,list->cols,ct->dim,ct->dims));
      for(j=1; j<list->rows; j++) {
//printf("%s - %s\n",list->data[i].desc,list->data[j].desc);
        buf[1] = list->data[j].val;
        fprintf(fo,",%.15lf",(sim->count==0)?func(buf,2,list->cols,ct->dim,ct->dims):1.0-func(buf,2,list->cols,ct->dim,ct->dims));
      }

      fprintf(fo,"\n");
    }

    ezgdal_show_progress(stdout,100,100);


    fclose(fo);
    sml_free_cell_type(ct);
    free_data(list);
    fclose(f);

  return 0;

}

