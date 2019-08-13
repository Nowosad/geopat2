#ifndef _LOCAL_PROTO_H_
#define _LOCAL_PROTO_H_

/****************************************************************************
 *
 * MODULE:	segmentation of motifels grid
 * AUTHOR(S):	Jaroslaw Jasiewicz, Jacek Niesterowicz, Tomasz Stepinski
 * PURPOSE:	information retrieval using categorical maps:
 *		compares grid of histograms
 * COPYRIGHT:	(C) Space Informatics Lab, University of Cincinnati
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

typedef enum {OUTPUT_LINKAGE,OPTION_REDUCTION,OPTION_LTHRESHOLD,OPTION_UTHRESHOLD,OPTION_SWAP,OPTION_WEIGHTS,OPTION_MEASURES,NUM_OF_LOCAL_OPTIONS} LOCAL_OPTIONS;
typedef enum {FLAG_SMOOTHING,FLAG_COMPLETE,FLAG_GROWING,FLAG_HIERARCHICAL,FLAG_THRESHOLD,FLAG_ALL,FLAG_QUAD,NUM_OF_LOCAL_FLAGS} LOCAL_FLAGS;
typedef enum {NORTH,SOUTH,WEST,EAST} DIRECTIONS;

struct area {
	/* part of the data which changes during segmantation process */
	double* histogram; /* contains centroid */
	int num_of_areas;
	struct fifo* hex_ids;
	int modified;
	double similarity_threshold; /* local homogenity */
	double homogenity; /* acually heterogeneity */
	double isolation;
	int id;
	int num_of_local_neighbors;
	int* local_neighborhood;

	/* ========================= */
	struct fifo* neighbors;
	int num_of_neighbors;
	int done;

	/* =====Geometry - for future use/remove === */
	double edge;
	int perimeter; /* in areas */
	int* bbox; /* in cells */
	int bbox_area;
	int bbox_perimeter;
};

struct diagnostics {
	double* cells_disagreement;
	double* segment_heterogeneity;
	double* segment_isolation;
	double* segment_threshold;
	int* index;
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
	/* recalculate distances to quantiles */
	double** tangents;
	double** quantiles;
	
} HEXGRID;

typedef struct {
	int minarea;
	int radius; /* size of the radius for homogenity */
	int window; /* size of the window for homogenity */
	int reduction;
	int use_smoothing;
	int complete_linkage;
	int all_layers; /* multilayer mode use _all_layers instead of _weighted_average */
	int quad_mode; /* type of topology */
	double null_threshold;
	double lower_similarity_threshold;
	double upper_similarity_threshold;
	double swap_threshold;
	int sampling_threshold; /* max number of elements in the matrix to start sampling */
	distance_func* calculate;
	S_PARAMS* parameters;
} LOCAL_PARAMS;

typedef struct {
	int index;
	int result; /* segment */ //TEN FRAGEMENT musi zostać ujednolicony: source, target, result segment itp.
	int target;
	double distance;
	double homogenity;
	int length_of_edge;
	int modified;
} DIST;

/* calculate */
int compare_by_dist (const void *a, const void *b);
int compare_by_dist_desc (const void *a, const void *b);
int sort_desc (const void * a, const void * b);
int sort_asc (const void * b, const void * a);
int sort_double_desc (const void * a, const void * b);
int sort_double_asc (const void * a, const void * b);
int* init_results(DATAINFO* d, LOCAL_PARAMS* p);
struct area** init_areas(DATAINFO* d, LOCAL_PARAMS* p);
double* new_histogram(DATAINFO* d, LOCAL_PARAMS* pars, unsigned index);
struct area* new_area(DATAINFO* d, LOCAL_PARAMS* pars, unsigned index);
int remove_area(struct area** a);
double* use_histogram(DATAINFO* d, int index);
double calculate_similarity(DATAINFO* d, LOCAL_PARAMS* p, struct area* a, struct area* b);
int read_histograms_to_memory(DATAINFO* d, LOCAL_PARAMS* p);
int* sample_histogram_ids(struct fifo* queue , int num_of_samples);
int add_histograms(DATAINFO* d, double* o, double* h, int num_of_areas);
int compare_grids_datainfo(DATAINFO* a, DATAINFO* b);
double interpolate(HEXGRID* hx, int layer, double distance);
double calculate2(HEXGRID* hx, LOCAL_PARAMS* p, double** pair);
int get_num_of_grids(char* files[]);

/* seeds */
int find_fisher_cut(DIST** distances,int num_of_dists,int first);
double* create_thresholds_map(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas);
int hex_find_quantiles(HEXGRID* hx, LOCAL_PARAMS* p, struct area** areas);

/* menu */
/*
struct Option* add_local_menu_item(LOCAL_OPTIONS option_name);
struct Flag* add_local_menu_flag(LOCAL_FLAGS flag_name);
*/

/* homogenity */
double** recalculate_histograms(DATAINFO* d, LOCAL_PARAMS* p, struct area** areas);
double* create_homogenity_map(DATAINFO* d, LOCAL_PARAMS* p, struct area** areas);
int find_edges(DATAINFO* d, LOCAL_PARAMS* p, struct area** areas);
int rebuild_histograms(DATAINFO* d, LOCAL_PARAMS* p, struct area** areas, double** orginal_histograms);
int build_init_histograms(DATAINFO* d, LOCAL_PARAMS* p, struct area** areas, double** orginal_histograms);

/*swap */
int swap_areas(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results);
double hex_calculate_linkage(HEXGRID* hx, LOCAL_PARAMS* p, struct area* a, struct area* b);
int hex_refine(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* p);

/* main functions */
unsigned* hex_find_seeds(HEXGRID* hx, LOCAL_PARAMS* p, struct area** areas, unsigned* num_of_seeds);
int hex_region_growing(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results, unsigned* seeds, unsigned num_of_seeds);
int hex_hierarchical(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results);
int hex_minarea (HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results);

/* in join */
int hex_join_areas(struct area* a, struct area* b);

/* write */
//int write_map(DATAINFO* d, LOCAL_PARAMS* p, void* map, char* output_map_name, int map_type, int pallete);
void write_raster(SML_DATA_HEADER* dh, void *map, char *raster_fname, int size);
void convert_to_vector(char *raster_fname, char *vector_fname);
int write_histogram_file(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* p, struct area** areas, char* filename);
int write_seed_layer(DATAINFO* d, int* seeds, int nseeds, char* filename);

/* diagnostics */
int calculate_segment_heterogeneity(HEXGRID* hx, LOCAL_PARAMS* pars, struct area* a);
double* save_thresholds(HEXGRID* hx, LOCAL_PARAMS* p, struct area** areas);

/* hex */
HEXGRID* hex_build_topology(DATAINFO** d, LOCAL_PARAMS* p, int num_of_layers, char* weights);
struct area** hex_build_areas(DATAINFO** d, HEXGRID* hx, LOCAL_PARAMS* p);
int* hex_remove_hexgrid(HEXGRID* hx);
int* hex_init_results(HEXGRID* hx);
double* hex_use_histogram(HEXGRID* hx, int index);
int hex_build_neighbors_queue(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results, int remove);
int hex_recalculate_neighbors(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int* results);
int hex_join_areas(struct area* a, struct area* b);
struct area* hex_new_area(HEXGRID* hx, unsigned index);
int hex_reclass(HEXGRID* hx, struct area** areas);
int* hex_create_segment_map(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* p, struct area** areas);
int hex_rebuild_histograms(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas);

/* local menu */
double* parse_weights(int num_of_layers, char* weights);

#endif
