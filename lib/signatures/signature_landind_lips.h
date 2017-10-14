#ifndef __LIPS_H__
#define __LIPS_H__

/****************************************************************************
 *
 * MODULE:	landscape indices signature
 * AUTHOR(S):	Jacek Niesterowicz
 * PURPOSE:	calculating a vector of landscape indices:
 *
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
#include <math.h>
#include <ctype.h>
#include <ezgdal.h>

/* list of categories in the input map */
typedef struct {
	int num;
	int* cat;
} MAP_CATS;

/* ============================ PARAMETERS =========================================*/

typedef struct {
/* parameters passed to histogram */
	int irregular_regions;
	/* LI: additional parameters */
	int eight_flag; /* 8-neighbor connectivity - 1; 4-neighbor - 0 */
	char** li_class_names; /* user-specified input indices */
	int li_num_of_class; /* number of indices to be calculated */
	int* li_class_flags; /* calculate an index or not */
	char** li_land_names;
	int* li_land_flags;
	int li_num_of_land;
	int* classes; /* list of classes to be calculated */
	int num_of_classes; /* # of classes ^ */
} H_PARAMS;

/* =========================== LANDSCAPE INDICES ============================ */
/* word "clump" means the same as "patch" in FRAGSTATS */
/* word "category" means the same as "class" in FRAGSTATS */

/* level coding */
#define ALL 0
#define CLASS 1
#define LANDSCAPE 2

typedef struct {
	unsigned long int id;
	int category;
	unsigned long int cells;
	double area; /* in map units */
	double perimeter; /* in map units */
	double centroid_r,centroid_c; /* as row/col coordinates */
	double gyrate;
	double contig;
	/*double linear; */
	/*double circle; */
} LI_CLUMP;

typedef struct {
	double total_area; /* including nulls; in map units */
	double notnull_area; /* area covered by not-null cells; in map units */
	double* cat_areas; /* in map units */
	double total_edge; /* including null edges, but not NOT_IN_REGION; in map units */
	double* cat_edges; /* total category edge including null-edges; in map units */
	double* cat_edge_matrix; /* between-category edge matrix; in map units */
	unsigned long int* adjacency; /* category adjacency matrix */
	unsigned long int total_adj; /* including class-null and class-NOT_IN_REGION adj. */
	unsigned long int* cat_adj; /* total adjacencies of categories */
	unsigned long int total_clumps; /* number of not-null clumps */
	unsigned long int* cat_clumps; /* number of clumps of each category */
	unsigned int ncats;
} LI_LANDSCAPE;

/* declaration of function type for landscape indices */
typedef double calculate_index(EZGDAL_LAYER*, LI_CLUMP*, LI_LANDSCAPE*, long int, long int, double, MAP_CATS*, int*);

typedef struct {
	calculate_index *index;
	char *name;
	char *description;
	int class_level;
	int landscape_level;
} LI_MENU;

/* ========================================================================== */

/* indices.c (misc functions) */
char* indices_menu_list(int);
int li_indices_number(int);
calculate_index* get_index_func(int);
char* get_index_name(int);
void li_set_params_all(H_PARAMS*, EZGDAL_LAYER*, int); /* hardcoded all indices */
void li_set_params_short(H_PARAMS*, EZGDAL_LAYER*, int); /* hardcoded similarity indices */
void li_free_parameters(H_PARAMS*);
int li_init_landscape(EZGDAL_LAYER*, LI_LANDSCAPE*);
void li_free_landscape(LI_LANDSCAPE*);
LI_CLUMP* li_add_clump(LI_CLUMP*, LI_LANDSCAPE*, unsigned long int, int, long int, long int);
void li_update_clump(LI_CLUMP*, LI_LANDSCAPE*, unsigned long int, int, long int, long int, double, int);
int li_check_gyrate(H_PARAMS*);
int li_check_contig(H_PARAMS*);
void li_gyrate(unsigned long int*, LI_CLUMP*, unsigned long int, long int, long int, double);
void li_contig(unsigned long int*, LI_CLUMP*, unsigned long int, long int, long int, double);
void li_write_names_to_header(char*, H_PARAMS*);

#endif
