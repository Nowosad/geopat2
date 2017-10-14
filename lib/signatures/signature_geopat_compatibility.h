#ifndef _COMPATIBILITY_H_
#define _COMPATIBILITY_H_

/****************************************************************************
 *
 * MODULE:	Compatibility with GRASS GeoPAT
 * AUTHOR(S):	Pawel Netzel
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <stdio.h>
#include <sml.h>

#define CELL_TYPE 0

typedef int CELL;
typedef double DCELL;
typedef int RASTER_MAP_TYPE;


typedef struct {
	char name[100];
	int *num_of_features;
	int *num_of_nulls;
	int *histogram;
} HISTOGRAM;


typedef struct {
  int rows;
  int cols;
  int proj;
  int zone;
  double ew_res;
  double ns_res;
  double north;
  double south;
  double west;
  double east;
} Cell_head;


typedef struct {
  char data_name[200];
  char output_name[200];
  char method[20];
  char list_of_maps[3000];
  int size_of_histogram;
  int num_of_histograms;
  HISTOGRAM* histograms;
  int pattern_size; /* size of the window */
  double add;
  Cell_head cell_hd;
  int taken;
  FILE *fd;
  int* buffer;
  double** all_histograms;
  SML_DATA_HEADER *dh;
} DATAINFO;

typedef struct {
  int size_of_histogram;
  int return_measure; /* 1 - similarity, 0 - distance */
  /* for co-occurence/weighted JSD */
  int num_of_cats;
  /* for EMD function */
  double *dist_matrix;
  double max_dist;
  /* for DTW function */
  int vector_size;
} S_PARAMS;


#ifndef MAX
#define MAX(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define MIN(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })
#endif

#define _(ARG) (ARG)

void G_message(const char *msg, ...);
void G_fatal_error(const char *msg, ...);
void G_warning(const char *msg, ...);
void* G_malloc(size_t size);
void G_free(void* v);

int Rast_is_c_null_value(const CELL* cellVal);
void Rast_set_c_null_value(CELL* cellVals, int numVals);
void Rast_set_d_null_value(DCELL* dcellVals, int numVals);


int* create_random_sequence(int length_of_sequence, int max_number);
char* check_input_names(char* name, const char* ext);

double* parse_weights(int num_of_layers, char* weights);
S_PARAMS* init_measure_parameters(int size_of_histogram, int measure);
int init_grid_datainfo(DATAINFO* d, char* filename, char* outputname);
int read_signatures_to_memory(DATAINFO* d);

#endif