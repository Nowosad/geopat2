#ifndef _LOCAL_PROTO_H_
#define _LOCAL_PROTO_H_

/****************************************************************************
 *
 * MODULE:	segmentation quality of motifels grid
 * AUTHOR(S):	Jaroslaw Jasiewicz, Jacek Niesterowicz, Tomasz Stepinski
 * PURPOSE:	information retrival using categorical maps:
 *		calculates quality of segmentation
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/

#include <math.h>
#include <string.h>
#include "macro.h"
#include "compatibility.h"
#include "list.h"

#include "../../lib/argtable/argtable3.h"
#include "../../lib/measures/measures.h"
#include "../../lib/signatures/signatures.h"
#include "../../lib/tools/libtools.h"


/* ======================================================= */

struct area {
	/* contains segment's info */
	double* histogram; /* contains centroid */
	int num_of_areas;
	int* hex_ids;
	double heterogeneity; /* = inhomogeneity */
	double isolation;
	int id;
	int num_of_samples;
	int* sample_ids;
	int* neighbors;
	int num_of_neighbors;
	int* edges;
	int done;
};

typedef struct {
	int level; /* level of scale reduction for refinement only */
	int nrows;
	int ncols;
	int nareas;
	int reduction;
	int size_of_histogram; /* total size if multilayer */
	int size_of_supermotifel; /* supermotifel size num of small elements */
/* part of the hex topology which is unmutable during segmentation process */
	int* low_res_grid;  /* for refinement only */
	double* thresholds;
	int num_of_subhistograms;
	int* sh_size_of_histogram; /* lengths of subistograms */
	double* sh_weights; /* weights of subhistograms */
	int** hex_neigborhoods;  /* for brick topology only */
	int** histogram_ids;
	double** histograms;
} HEXGRID;

typedef struct {
	int reduction;
	int complete_linkage;
	int quad_mode; /* type of topology */
	int no_weight;
	double null_threshold;
	int sampling_threshold; /* max number of elements in the matrix to start sampling */
	distance_func* calculate;
	S_PARAMS* parameters;
} LOCAL_PARAMS;

/* calculate */
int sort_desc (const void * a, const void * b);
int sort_asc (const void * b, const void * a);
int sort_double_desc (const void * a, const void * b);
int sort_double_asc (const void * a, const void * b);
double* use_histogram(DATAINFO* d, int index);
int add_histograms(DATAINFO* d, double* o, double* h, int num_of_areas);

/* write */
void write_raster(SML_DATA_HEADER* dh, void *map, char *raster_fname, int size);
void write_raster_dbl(SML_DATA_HEADER* dh, void *map, char *raster_fname, int size);

/* diagnostics */
double* calculate_heterogeneity(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas);
double* calculate_isolation(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas);

/* rebuild */
struct area** build_areas(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* p, int* results);
int build_neighbors_queue(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results);
int* create_segment_map(HEXGRID* hx, LOCAL_PARAMS* p, EZGDAL_LAYER* l);

/* hex */
HEXGRID* hex_build_topology(DATAINFO* d, LOCAL_PARAMS* p);
int hex_remove_hexgrid(HEXGRID* hx);
double* hex_use_histogram(HEXGRID* hx, int index);
struct area* hex_new_area(HEXGRID* hx, unsigned index);
int remove_areas(struct area** areas, int num_of_areas);

#endif
