/****************************************************************************
 *
 * MODULE:	Compatibility with GRASS GeoPAT
 * AUTHOR(S):	Pawel Netzel, Jakub Nowosad
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <math.h>

#include "../../lib/ezGDAL/ezgdal.h"
#include "../../lib/SML/sml.h"

#include "compatibility.h"

#define WARN 0
#define ERR  1
#define MSG  2

double Max(double a, double b) {
	if(a>b)
		return a;
	else
		return b;
}

double Min(double a, double b) {
	if(a>b)
		return b;
	else
		return a;
}


static void print_msg(int type, const char *template, va_list ap) {
  char buffer[2000];

  if(type==ERR)
    sprintf(buffer, "ERROR: %s\n", template);
  else if(type==WARN)
    sprintf(buffer, "WARNING: %s\n", template);
  else
    sprintf(buffer, "%s\n", template);

  vprintf(buffer, ap);
}

void G_message(const char *msg, ...) {
  va_list ap;

  va_start(ap, msg);
  print_msg(MSG, msg, ap);
  va_end(ap);

}

void G_fatal_error(const char *msg, ...) {
  static int busy;
  va_list ap;

  if (busy)
    exit(1);
  busy = 1;

  va_start(ap, msg);
  print_msg(ERR, msg, ap);
  va_end(ap);

  exit(1);
}

void G_warning(const char *msg, ...) {
  va_list ap;

  va_start(ap, msg);
  print_msg(WARN, msg, ap);
  va_end(ap);
}

void* G_malloc(size_t size) {
  return malloc(size);
}

void G_free(void* v) {
  free(v);
}

int Rast_is_c_null_value(const CELL * cellVal) {
    return *cellVal == (CELL) 0x80000000;
}

void Rast_set_c_null_value(CELL * cellVals, int numVals) {
    int i;

    for (i = 0; i < numVals; i++)
	cellVals[i] = (int)0x80000000;
}

void Rast_set_d_null_value(DCELL * dcellVals, int numVals) {
    static const unsigned char null_bits[8] = {
	0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF
    };
    int i;

    for (i = 0; i < numVals; i++)
	memcpy(&dcellVals[i], null_bits, sizeof(null_bits));
}



/***********************************************/

int* create_random_sequence(int length_of_sequence, int max_number) {
	int i, j, tmp;
	int* randoms;
	int* ret;

	srand((int)time(NULL)); /* seed in the future */
	randoms=malloc(max_number*sizeof(int));
	for(i=0;i<max_number;++i)
		randoms[i]=i+1;

	/* shuffle */
	for(i=0;i<max_number;++i) {
		j = i + rand() / (RAND_MAX / (max_number - i) + 1);
		tmp=randoms[j];
		randoms[j]=randoms[i];
		randoms[i]=tmp;
	}

//	qsort(randoms, length_of_sequence, sizeof(int), compare_by_int);

	ret=malloc((length_of_sequence+1)*sizeof(int));
	for(i=0;i<length_of_sequence;++i)
		ret[i]=randoms[i]-1;
	ret[length_of_sequence]=0;
	free(randoms);
	return ret;
}

char* check_input_names(char* name, const char* ext) { return "";}

S_PARAMS* init_measure_parameters(int size_of_histogram, int measure) {

    S_PARAMS *p = malloc(sizeof(S_PARAMS));

    p->size_of_histogram=size_of_histogram;
    p->return_measure=measure;

    return p;
}



int init_grid_datainfo(DATAINFO* d, char* filename, char* outputname) {

    SML_DATA_HEADER  *dh = sml_open_layer(filename);

    if(!dh)
        G_fatal_error(_("Cannot open signatures file %s"), filename);

    d->dh = dh;

    strcpy(d->data_name,filename);
    strcpy(d->output_name,outputname);

    d->list_of_maps[0] = '\0';
    d->method[0] = '\0';
    d->size_of_histogram = dh->cell_N_elements;
    d->cell_hd.rows = dh->file_win->rows;
    d->cell_hd.cols = dh->file_win->cols;
    d->cell_hd.proj = 0;
    d->cell_hd.zone = 0;
    d->cell_hd.ew_res = fabs(dh->file_win->at[1]);
    d->cell_hd.ns_res = fabs(dh->file_win->at[5]);
    d->cell_hd.north = dh->file_win->at[3];
    d->cell_hd.south = d->cell_hd.north - d->cell_hd.rows * d->cell_hd.ns_res;
    d->cell_hd.west = dh->file_win->at[0];
    d->cell_hd.east = d->cell_hd.west + d->cell_hd.cols * d->cell_hd.ew_res;
    d->pattern_size = 1;
    d->buffer=NULL; /*malloc(nhists*(d->size_of_histogram+2)*sizeof(int));*/

        /* allocate histograms... */
    d->histograms = NULL;
    /*d->histograms=malloc(nhists*sizeof(HISTOGRAM));*/


    return 0;
}

double* parse_weights(int num_of_layers, char* weights)
{
    /* parser */
    
    if(num_of_layers==1 && weights) {
        G_warning("Ignore weights for a single-layer segmentation");
        return NULL;
    }
    
    if(num_of_layers==1)
        return NULL;
    
    int i;
    double* parsed=calloc(num_of_layers,sizeof(double));
    if(num_of_layers>1 && !weights) { /* set weights to equal 1/n */
    for(i=0;i<num_of_layers;++i)
        parsed[i]=1./(double)num_of_layers;
        return parsed;
    }
    
    char buf[10000];
    strcpy(buf,weights);
    char delim[2]=",";
    char* token;
    
    token=strtok(buf,delim);
    int col=0;
    
    while(token) {
        if(col>num_of_layers)
            G_fatal_error("Too many weights. Expected: %d",num_of_layers);
        parsed[col]=atoi(token);
        token=strtok(NULL,delim);
        col++;
    }
    if(col>num_of_layers)
        G_fatal_error("Too few weights. Expected: %d",num_of_layers);
    
    /* recalculate so we can put any list of numbers to estabilish proportions */
    double sum=0;
    
    for(i=0;i<num_of_layers;++i) {
        if(parsed[i]<=0)
            G_fatal_error("Sum of weights cannot be zero");
        sum+=parsed[i];
    }
    
    if(sum==0)
        G_fatal_error("Sum of weights cannot be zero");
    
    for(i=0;i<num_of_layers;++i)
        parsed[i]/=sum;
    
    return parsed;
}


int read_signatures_to_memory(DATAINFO* d) {
    int r,c, size;
    int nrows=d->cell_hd.rows;
    int ncols=d->cell_hd.cols;
    long int ncells=nrows*ncols;
    long int i = 0;
    void *row, *cell;

/*    printf("Reading data: ...    0%%");*/
    ezgdal_show_progress(stdout,0,nrows);
    d->all_histograms=malloc(ncells*sizeof(double*));
    size = d->size_of_histogram*sizeof(double);
    row = sml_create_cell_row_buffer(d->dh);

    for(r=0; r<nrows; r++) {

        ezgdal_show_progress(stdout,r,nrows);
        sml_read_row_from_layer(d->dh,row,r);
        for(c=0; c<ncols; c++) {
            cell = sml_get_cell_pointer(d->dh, row, c);
            if(sml_is_cell_null(cell))
                d->all_histograms[i] = NULL;
            else {
                d->all_histograms[i] = malloc(size);
                memcpy(d->all_histograms[i], sml_get_cell_data(cell), size);
            }
            i++;
        }

    }
    ezgdal_show_progress(stdout,nrows,nrows);

    free(row);

    return 0;

}

