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


/******************************************************************************
Landscape indices (metrics) functions.
Names are the shortcuts used in the FRAGSTATS software
word "clump" means the same as "patch" in FRAGSTATS 
word "category" means the same as "class" in FRAGSTATS 
Author: Jacek Niesterowicz, University of Cincinnati
******************************************************************************/

#include <float.h>
#include "signature_landind_lips.h"
#include "signature_landind_macro.h"
#include "signature_geopat_compatibility.h"

extern calculate_index pland;
extern calculate_index lpi;
extern calculate_index ed;
extern calculate_index area_mn;
extern calculate_index area_am;
extern calculate_index area_md;
extern calculate_index area_ra;
extern calculate_index area_sd;
extern calculate_index area_cv;
extern calculate_index gyrate_mn;
extern calculate_index gyrate_am;
extern calculate_index gyrate_md;
extern calculate_index gyrate_ra;
extern calculate_index gyrate_sd;
extern calculate_index gyrate_cv;
extern calculate_index pafrac;
extern calculate_index para_mn;
extern calculate_index para_am;
extern calculate_index para_md;
extern calculate_index para_ra;
extern calculate_index para_sd;
extern calculate_index para_cv;
extern calculate_index shape_mn;
extern calculate_index shape_am;
extern calculate_index shape_md;
extern calculate_index shape_ra;
extern calculate_index shape_sd;
extern calculate_index shape_cv;
extern calculate_index frac_mn;
extern calculate_index frac_am;
extern calculate_index frac_md;
extern calculate_index frac_ra;
extern calculate_index frac_sd;
extern calculate_index frac_cv;
/*
extern calculate_index linear_mn;
extern calculate_index linear_am;
extern calculate_index linear_md;
extern calculate_index linear_ra;
extern calculate_index linear_sd;
extern calculate_index linear_cv;
extern calculate_index circle_mn;
extern calculate_index circle_am;
extern calculate_index circle_md;
extern calculate_index circle_ra;
extern calculate_index circle_sd;
extern calculate_index circle_cv;
*/
extern calculate_index contig_mn;
extern calculate_index contig_am;
extern calculate_index contig_md;
extern calculate_index contig_ra;
extern calculate_index contig_sd;
extern calculate_index contig_cv;
extern calculate_index contag;
extern calculate_index iji;
extern calculate_index pladj;
extern calculate_index ai;
extern calculate_index lsi;
extern calculate_index cohesion;
extern calculate_index pd;
extern calculate_index split;
extern calculate_index division;
extern calculate_index mesh;
extern calculate_index prd;
extern calculate_index rpr;
extern calculate_index shdi;
extern calculate_index sidi;
extern calculate_index msidi;
extern calculate_index shei;
extern calculate_index siei;
extern calculate_index msiei;

typedef enum {
li_pland,
li_lpi,
li_ed,
li_area_mn,
li_area_am,
li_area_md,
li_area_ra,
li_area_sd,
li_area_cv,
li_gyrate_mn,
li_gyrate_am,
li_gyrate_md,
li_gyrate_ra,
li_gyrate_sd,
li_gyrate_cv,
li_pafrac,
li_para_mn,
li_para_am,
li_para_md,
li_para_ra,
li_para_sd,
li_para_cv,
li_shape_mn,
li_shape_am,
li_shape_md,
li_shape_ra,
li_shape_sd,
li_shape_cv,
li_frac_mn,
li_frac_am,
li_frac_md,
li_frac_ra,
li_frac_sd,
li_frac_cv,
/*
li_linear_mn,
li_linear_am,
li_linear_md,
li_linear_ra,
li_linear_sd,
li_linear_cv,
li_circle_mn,
li_circle_am,
li_circle_md,
li_circle_ra,
li_circle_sd,
li_circle_cv,
*/
li_contig_mn,
li_contig_am,
li_contig_md,
li_contig_ra,
li_contig_sd,
li_contig_cv,
li_contag,
li_iji,
li_pladj,
li_ai,
li_lsi,
li_cohesion,
li_pd,
li_split,
li_division,
li_mesh,
li_prd,
li_rpr,
li_shdi,
li_sidi,
li_msidi,
li_shei,
li_siei,
li_msiei,
num_of_indices } INDICES;

static LI_MENU li_menu[] = {

	/* simple composition of categories */
	{pland,"pland","class percentage of landscape",1,0},

	/* Area/Edge metrics */
	{lpi,"lpi","largest patch index",1,1},
	{ed,"ed","edge density",1,1},
	{area_mn,"area_mn","patch area mean",1,1},
	{area_am,"area_am","patch area area-weighted mean",1,1},
	{area_md,"area_md","patch area median",1,1},
	{area_ra,"area_ra","patch area range",1,1},
	{area_sd,"area_sd","patch area standard deviation",1,1},
	{area_cv,"area_cv","patch area coeff. of variation",1,1},
	{gyrate_mn,"gyrate_mn","radius of gyration mean",1,1},
	{gyrate_am,"gyrate_am","radius of gyration area-weighted mean",1,1},
	{gyrate_md,"gyrate_md","radius of gyration median",1,1},
	{gyrate_ra,"gyrate_ra","radius of gyration range",1,1},
	{gyrate_sd,"gyrate_sd","radius of gyration standard deviation",1,1},
	{gyrate_cv,"gyrate_cv","radius of gyration coeff. of variation",1,1},

	/* Shape metrics */
	{pafrac,"pafrac","perimeter-area fractal dimension",1,1},
	{para_mn,"para_mn","perimeter-area ratio mean",1,1},
	{para_am,"para_am","perimeter-area ratio area-weighted mean",1,1},
	{para_md,"para_md","perimeter-area ratio median",1,1},
	{para_ra,"para_ra","perimeter-area ratio range",1,1},
	{para_sd,"para_sd","perimeter-area ratio standard deviation",1,1},
	{para_cv,"para_cv","perimeter-area ratio coeff. of variation",1,1},
	{shape_mn,"shape_mn","shape index mean",1,1},
	{shape_am,"shape_am","shape index area-weighted mean",1,1},
	{shape_md,"shape_md","shape index median",1,1},
	{shape_ra,"shape_ra","shape index range",1,1},
	{shape_sd,"shape_sd","shape index standard deviation",1,1},
	{shape_cv,"shape_cv","shape index coeff. of variation",1,1},
	{frac_mn,"frac_mn","fractal index mean",1,1},
	{frac_am,"frac_am","fractal index area-weighted mean",1,1},
	{frac_md,"frac_md","fractal index median",1,1},
	{frac_ra,"frac_ra","fractal index range",1,1},
	{frac_sd,"frac_sd","fractal index standard deviation",1,1},
	{frac_cv,"frac_cv","fractal index coeff. of variation",1,1},
/*
	{linear_mn,"linear_mn","linearity index mean"},
	{linear_am,"linear_am","linearity index area-weighted mean"},
	{linear_md,"linear_md","linearity index median"},
	{linear_ra,"linear_ra","linearity index range"},
	{linear_sd,"linear_sd","linearity index standard deviation"},
	{linear_cv,"linear_cv","linearity index coeff. of variation"},
	{circle_mn,"circle_mn","related circumscribing circle mean"},
	{circle_am,"circle_am","related circumscribing circle area-weighted mean"},
	{circle_md,"circle_md","related circumscribing circle median"},
	{circle_ra,"circle_ra","related circumscribing circle range"},
	{circle_sd,"circle_sd","related circumscribing circle standard deviation"},
	{circle_cv,"circle_cv","related circumscribing circle coeff. of variation"},
*/
	{contig_mn,"contig_mn","contiguity index mean",1,1},
	{contig_am,"contig_am","contiguity index area-weighted mean",1,1},
	{contig_md,"contig_md","contiguity index median",1,1},
	{contig_ra,"contig_ra","contiguity index range",1,1},
	{contig_sd,"contig_sd","contiguity index standard deviation",1,1},
	{contig_cv,"contig_cv","contiguity index coeff. of variation",1,1},

	/* Aggregation metrics */
	{contag,"contag","contagion index",0,1},
	{iji,"iji","interspersion & juxtaposition index",1,1},
	{pladj,"pladj","percentage of like adjacencies",1,1},
	{ai,"ai","aggregation index",1,1},
	{lsi,"lsi","landscape shape index",1,1},
	{cohesion,"cohesion","patch cohesion index",1,1},
	{pd,"pd","patch density",1,1},
	{split,"split","splitting index",1,1},
	{division,"division","landscape division index",1,1},
	{mesh,"mesh","effective mesh size",1,1},
	
	/* Diversity metrics */
	{prd,"prd","patch richness density",0,1},
	{rpr,"rpr","relative patch richness",0,1},
	{shdi,"shdi","Shannon's diversity index",0,1},
	{sidi,"sidi","Simpson's diversity index",0,1},
	{msidi,"msidi","modified Simpson's diversity index",0,1},
	{shei,"shei","Shannon's evenness index",0,1},
	{siei,"siei","Simpson's evenness index",0,1},
	{msiei,"msiei","modified Simpson's evenness index",0,1},

	{NULL,0,0,0,0}
};

/* ================ misc functions ===========================================*/

int li_indices_number(int level)
/*
level=0 - number of all indices
level=1 - number of class level indices
level=2 - number of landscape level indices
*/
{
	int i,count=0;

	switch(level){
		case ALL:
			return num_of_indices;
		case CLASS:
			for(i=0;i<num_of_indices;++i){
				if(li_menu[i].class_level)
					count++;
			}
			return count;
		case LANDSCAPE:
			for(i=0;i<num_of_indices;++i){
				if(li_menu[i].landscape_level)
					count++;
			}
			return count;
		default:
			G_fatal_error("indices_number() says: Wrong level id!");
	}
	return 0;
}

char* indices_menu_list(int level)
/*
level=0 - print all indices names;
level=1 - print class level indices names; 
level=2 - print landscape level indices names.
*/
{
	int i,count=0;
	int n=li_indices_number(level);
	static char* full_list;
	full_list=malloc(2000);
	for(i=0;li_menu[i].name;++i) {
		if(level==CLASS){
			if(!li_menu[i].class_level){
				continue;
			}
		}
		else if(level==LANDSCAPE){
			if(!li_menu[i].landscape_level){
				continue;
			}
		}
		count++;
		strcat(full_list,li_menu[i].name);
		strcat(full_list,(count<n)?",":"");
	}
	return full_list;
}

calculate_index* get_index_func(int number)
{
	return li_menu[number].index;
}

char* get_index_name(int number)
{
	return li_menu[number].name;
}

void li_set_params_all(H_PARAMS* p, EZGDAL_LAYER* l, int eight_flag){
	/* eight_flag: 0 - 4-neighborhood; 1 - 8-neighborhood */
	int count,i,j;

	p->eight_flag = eight_flag;
	/* pick all class level */
	p->li_class_flags=calloc(num_of_indices,sizeof(int)); /* init to 0 */
	count=0;
	p->li_num_of_class=li_indices_number(CLASS);
	p->li_class_names=malloc(p->li_num_of_class*sizeof(char*));
	for(i=0; i<num_of_indices; ++i)
		if(li_menu[i].class_level){
			p->li_class_flags[i]=1;
			p->li_class_names[count++]=li_menu[i].name;
		}
	p->num_of_classes=l->stats->map_max_val+1;
	p->classes=malloc(p->num_of_classes*sizeof(int));
	j=0;
	for(i=0; i<l->stats->hist_N; i++)
		if(l->stats->map_cat[i] >= 0)
			p->classes[j++] = i + (int)(l->stats->min);

	/* pick all landscape level */
	p->li_land_flags=calloc(num_of_indices,sizeof(int)); // init to 0
	count=0;
	p->li_num_of_land=li_indices_number(LANDSCAPE);
	p->li_land_names=malloc(p->li_num_of_land*sizeof(char*));
	for(i=0; i<num_of_indices; ++i)
		if(li_menu[i].landscape_level){
			p->li_land_flags[i]=1;
			p->li_land_names[count++]=li_menu[i].name;
		}
	return;
}

/* sets parameters for computing only composition and all landscape level indices */
void li_set_params_short(H_PARAMS* p, EZGDAL_LAYER* l, int eight_flag){
	/* eight_flag: 0 - 4-neighborhood; 1 - 8-neighborhood */
	int count,i,j;

	p->eight_flag = eight_flag;
	/* pick all class level */
	p->li_class_flags=calloc(num_of_indices,sizeof(int)); /* init to 0 */
	p->li_num_of_class=1;
	p->li_class_names=malloc(p->li_num_of_class*sizeof(char*));
	for(i=0; i<num_of_indices; ++i){
		if((strcmp(li_menu[i].name,"pland"))==0){
			p->li_class_flags[i]=1;
			p->li_class_names[0]="pland";
		}
	}
	p->num_of_classes=l->stats->map_max_val+1;
	p->classes=malloc(p->num_of_classes*sizeof(int));
	j=0;
	for(i=0; i<l->stats->hist_N; i++)
		if(l->stats->map_cat[i] >= 0)
			p->classes[j++] = i + (int)(l->stats->min);

	/* pick all landscape level */
	p->li_land_flags=calloc(num_of_indices,sizeof(int)); // init to 0
	count=0;
	p->li_num_of_land=li_indices_number(LANDSCAPE);
	p->li_land_names=malloc(p->li_num_of_land*sizeof(char*));
	for(i=0; i<num_of_indices; ++i)
		if(li_menu[i].landscape_level){
			p->li_land_flags[i]=1;
			p->li_land_names[count++]=li_menu[i].name;
		}
	return;
}

void li_free_parameters(H_PARAMS* p)
{
	if(p->li_class_names)
		free(p->li_class_names);
	if(p->li_land_names)
		free(p->li_land_names);
	if(p->li_class_flags)
		free(p->li_class_flags);
	if(p->li_land_flags)
		free(p->li_land_flags);
	if(p->classes)
		free(p->classes);
	return;
}

int li_init_landscape(EZGDAL_LAYER* l, LI_LANDSCAPE* landscape)
{
	int length = (int)l->stats->max + 1; // highest category number + 1

	/* initialize landscape statistics variables */
	landscape->total_area = 0;
	landscape->notnull_area = 0;
	landscape->cat_areas = calloc(length, sizeof(double));
	landscape->total_edge = 0;
	landscape->cat_edges = calloc(length, sizeof(double));
	landscape->cat_edge_matrix = calloc(length*length, sizeof(double));
	landscape->adjacency = calloc(length*length, sizeof(unsigned long int));
	landscape->total_adj = 0;
	landscape->cat_adj = calloc(length, sizeof(unsigned long int));
	landscape->total_clumps = 0;
	landscape->cat_clumps = calloc(length, sizeof(unsigned long int));
	landscape->ncats = 0;
	return length;
}

void li_free_landscape(LI_LANDSCAPE* landscape)
{
	/* free allocated arrays */
	free(landscape->cat_areas);
	free(landscape->cat_edges);
	free(landscape->cat_edge_matrix);
	free(landscape->adjacency);
	free(landscape->cat_adj);
	free(landscape->cat_clumps);
	return;
}

LI_CLUMP* li_add_clump(LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	unsigned long int id, int cat, long int row, long int col)
{
	/* function runs when a new clump is found */

	/* initialize clump statistics */
	cl_stats = realloc(cl_stats, (id+1)*sizeof(LI_CLUMP));
	cl_stats[id].id = id;
	cl_stats[id].category = cat;
	cl_stats[id].cells = 0;
	cl_stats[id].area = 0; 
	cl_stats[id].perimeter = 0;
	cl_stats[id].centroid_r = row;
	cl_stats[id].centroid_c = col;
	cl_stats[id].gyrate = 0;
	cl_stats[id].contig = 0;
	//cl_stats[id].linear = 0;
	//cl_stats[id].circle = 0;

	/* update landscape statistics */
	if(!(la_stats->cat_clumps[cat]))
		la_stats->ncats++;
	la_stats->total_clumps++;
	la_stats->cat_clumps[cat]++;
	
	return cl_stats;
}

void li_update_clump(LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats, 
	unsigned long int id, int cat, long int row, long int col, double res_area, int eight_flag)
{
	/* update area statistics */
	cl_stats[id].cells++;
	cl_stats[id].area+=res_area;
	/* recalc centroid */
	cl_stats[id].centroid_r = cl_stats[id].centroid_r + ((double)row-cl_stats[id].centroid_r)/cl_stats[id].cells;
	cl_stats[id].centroid_c = cl_stats[id].centroid_c + ((double)col-cl_stats[id].centroid_c)/cl_stats[id].cells;
	/* update landscape statistics */
	la_stats->notnull_area+=res_area;
	la_stats->cat_areas[cat]+=res_area;
	la_stats->total_adj+=4; // 4 or 8 adjacent cells of current not-null cell
	la_stats->cat_adj[cat]+=4;
	return;
}

/* check if radius of gyration is to be calculated */
int li_check_gyrate(H_PARAMS* p)
{
	if(p->li_land_flags[li_gyrate_mn]) return 1;
	if(p->li_land_flags[li_gyrate_am]) return 1;
	if(p->li_land_flags[li_gyrate_md]) return 1;
	if(p->li_land_flags[li_gyrate_ra]) return 1;
	if(p->li_land_flags[li_gyrate_sd]) return 1;
	if(p->li_land_flags[li_gyrate_cv]) return 1;
	if(p->li_class_flags[li_gyrate_mn]) return 1;
	if(p->li_class_flags[li_gyrate_am]) return 1;
	if(p->li_class_flags[li_gyrate_md]) return 1;
	if(p->li_class_flags[li_gyrate_ra]) return 1;
	if(p->li_class_flags[li_gyrate_sd]) return 1;
	if(p->li_class_flags[li_gyrate_cv]) return 1;
	return 0;
}

/* function calculates radius of gyration for each clump/patch */
void li_gyrate(unsigned long int* map_clump, LI_CLUMP* cl_stats, unsigned long int nclumps, long int nrows, long int ncols, double resolution)
{
	long int row,col;
	unsigned long int id;

	/* for each pixel in the landscape */
	for(row=0; row<nrows; ++row){
		for(col=0; col<ncols; ++col){
			if((id=map_clump[row*ncols+col])==0)
				continue;
			cl_stats[id].gyrate += sqrt((row-round(cl_stats[id].centroid_r))*(row-round(cl_stats[id].centroid_r))+
				(col-round(cl_stats[id].centroid_c))*(col-round(cl_stats[id].centroid_c)));
		}
	}
	/* mean distance */
	for(id=1; id<nclumps+1; ++id)
		cl_stats[id].gyrate = cl_stats[id].gyrate / cl_stats[id].cells * resolution;
	return;
}

/* check if contiguity index is to be calculated */
int li_check_contig(H_PARAMS* p)
{
	if(p->li_land_flags[li_contig_mn]) return 1;
	if(p->li_land_flags[li_contig_am]) return 1;
	if(p->li_land_flags[li_contig_md]) return 1;
	if(p->li_land_flags[li_contig_ra]) return 1;
	if(p->li_land_flags[li_contig_sd]) return 1;
	if(p->li_land_flags[li_contig_cv]) return 1;
	if(p->li_class_flags[li_contig_mn]) return 1;
	if(p->li_class_flags[li_contig_am]) return 1;
	if(p->li_class_flags[li_contig_md]) return 1;
	if(p->li_class_flags[li_contig_ra]) return 1;
	if(p->li_class_flags[li_contig_sd]) return 1;
	if(p->li_class_flags[li_contig_cv]) return 1;
	return 0;
}

/* contugity index for each clump */
void li_contig(unsigned long int* map_clump, LI_CLUMP* cl_stats, unsigned long int nclumps, long int nrows, long int ncols, double resolution)
{
	unsigned long int i,id;
	long int r,c,next_r,next_c;
	long int* sum = calloc((nclumps+1),sizeof(long int)); // sum of products for each clump
	int template[9] = { 1, 1, 2, 1, 2, 1, 2, 1, 2 }; // 3x3 kernel template
	int temp_sum=0;
	for(i=0; i<9; ++i)
		temp_sum+=template[i];

	/* for each pixel in the landscape */
	for(r=0; r<nrows; ++r){
		for(c=0; c<ncols; ++c){
			if((id=map_clump[r*ncols+c])==0)
				continue; // null

			/* check 9-cell neighborhood */
			/* template of neighborhood in macro.h */
			for(i=0; i<9; ++i){
				if(NOT_IN_REGION(i))
					continue;
				next_r = SNR(i);
				next_c = SNC(i);

				/* instead of binary map: if ids agree, add template value */
				if(id==map_clump[next_r*ncols+next_c])
					sum[id] += template[i];
			}
		}
	} // end row
	
	for(id=1; id<=nclumps; ++id)
		cl_stats[id].contig = ((double)sum[id]/cl_stats[id].cells-1.) / (double)(temp_sum-1.);
	free(sum);
	return;
}

void li_write_names_to_header(char* opt_output_name, H_PARAMS* p)
{
	FILE *header_fp;
	char header_name[500];
	int i;
	
	sprintf(header_name,"%s.hr",opt_output_name);
	header_fp=fopen(header_name, "a"); // open for appending
	fprintf(header_fp,"indices:");
	for(i=0; i<p->li_num_of_land; ++i)
		fprintf(header_fp,"%s%s",p->li_land_names[i],(i==(p->li_num_of_land-1))?"\n":",");
	fclose(header_fp);
	return;
}

/*========================LANDSCAPE INDICES FUNCTIONS=========================*/

/* internal function for qsort() */
int compare(const void* a, const void* b)
{
	return (int)((*(double*)a - *(double*)b) * 1000);
}

/*========================Area/Edge indices===================================*/

/* class percent of landscape (%) */
double pland(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	/* "pland" option calculates the metric for every category in the map */
	/* number of pland indices = number of categories in the map */
	if(!class){ // landscape level
		G_fatal_error(_("\"pland\" is class level only."));
	}
	else{ // single class level
		return la_stats->cat_areas[*class] / la_stats->notnull_area * 100;
	}
	return 0;
}

/* largest patch index (%) */
double lpi(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double max=0;
	
	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].area>max)
				max = cl_stats[i].area;
		}
	}
	else{ // single class level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class && cl_stats[i].area>max)
				max = cl_stats[i].area;
		}
	}
	return max / la_stats->notnull_area * 100;
}

/* edge density [m/ha] */
double ed(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	/* FRAGSTATS counts edges twice */
	
	if(!class){ // landscape level
		return la_stats->total_edge/la_stats->notnull_area*10000;
	}
	else{ // single class level
		return la_stats->cat_edges[*class]/la_stats->notnull_area*10000;
	}
}

/* patch area [ha] */
/* patch area mean */
double area_mn(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	
	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].area;
		return (sum / la_stats->total_clumps) / 10000;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here

		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].area;
		return (sum / la_stats->cat_clumps[*class]) / 10000;
	}
}

/* patch area area-weighted mean */
double area_am(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i; 
	double sum=0;
	
	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += (cl_stats[i].area * (cl_stats[i].area / la_stats->notnull_area)) / 10000;
		return sum;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += (cl_stats[i].area * (cl_stats[i].area / la_stats->cat_areas[*class])) / 10000;
		return sum;
	}
}

/* patch area median */
double area_md(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double median;
	double* list;
	
	if(!class){ // landscape level
		/* create temporary table and quicksort it */
		list = malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			list[i-1] = cl_stats[i].area;
		qsort(list, la_stats->total_clumps, sizeof(double), compare);

		if(la_stats->total_clumps%2)
			median = list[la_stats->total_clumps/2];
		else
			median = (list[la_stats->total_clumps/2]+list[la_stats->total_clumps/2-1])/2;
		free(list);
		return median / 10000;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		/* create temporary table and quicksort it */
		long int current=0;
		list = malloc(la_stats->cat_clumps[*class]*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				list[current++] = cl_stats[i].area;
		qsort(list, la_stats->cat_clumps[*class], sizeof(double), compare);

		if(la_stats->cat_clumps[*class]%2)
			median = list[la_stats->cat_clumps[*class]/2];
		else
			median = (list[la_stats->cat_clumps[*class]/2]+list[la_stats->cat_clumps[*class]/2-1])/2;
		free(list);
		return median / 10000;
	}
}

/* patch area range */
double area_ra(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double min=DBL_MAX;
	double max=0;

	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].area>max)
				max = cl_stats[i].area;
			if(cl_stats[i].area<min)
				min = cl_stats[i].area;
		}
		return (max-min) / 10000;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category!=*class)
				continue; // not this class
			if(cl_stats[i].area>max)
				max = cl_stats[i].area;
			if(cl_stats[i].area<min)
				min = cl_stats[i].area;
		}
		return (max-min) / 10000;
	}
}

/* patch area standard deviation */
double area_sd(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;

	if(!class){ // landscape level
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].area;
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			dev = cl_stats[i].area - mean;
			sd += dev*dev;
		}
		return sqrt(sd / la_stats->total_clumps) / 10000;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].area;
		mean = sum / la_stats->cat_clumps[*class];

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class){
				dev = cl_stats[i].area - mean;
				sd += dev*dev;
			}
		}
		return sqrt(sd / la_stats->cat_clumps[*class]) / 10000;
	}
}

/* patch area coefficient of variation */
double area_cv(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;

	if(!class){ // landscape level
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].area;
		mean = (double)sum / la_stats->total_clumps;

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			dev = cl_stats[i].area - mean;
			sd += dev*dev;
		}
		sd = sqrt(sd / la_stats->total_clumps);
		return sd / mean * 100;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].area;
		mean = sum / la_stats->cat_clumps[*class];

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class){
				dev = cl_stats[i].area - mean;
				sd += dev*dev;
			}
		}
		sd = sqrt(sd / la_stats->cat_clumps[*class]);
		return sd / mean * 100;	
	}
}

/* radius of gyration */
/* radius of gyration mean */
double gyrate_mn(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].gyrate;
		return sum / la_stats->total_clumps;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].gyrate;
		return (sum / la_stats->cat_clumps[*class]) / 10000;
	}
}

/* radius of gyration area-weighted mean */
double gyrate_am(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i; 
	double sum=0;
	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].gyrate * (cl_stats[i].area / la_stats->notnull_area);
		return sum;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += (cl_stats[i].gyrate * (cl_stats[i].area / la_stats->cat_areas[*class])) / 10000;
		return sum;
	}
}

/* radius of gyration median */
double gyrate_md(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double median;
	double* list;

	if(!class){ // landscape level
		/* create temporary table and quicksort it */
		list = malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			list[i-1] = cl_stats[i].gyrate;
		qsort(list, la_stats->total_clumps, sizeof(double), compare);

		if(la_stats->total_clumps%2)
			median = list[la_stats->total_clumps/2];
		else
			median = (list[la_stats->total_clumps/2]+list[la_stats->total_clumps/2-1])/2;
		free(list);
		return median;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		/* create temporary table and quicksort it */
		long int current=0;
		list = malloc(la_stats->cat_clumps[*class]*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				list[current++] = cl_stats[i].gyrate;
		qsort(list, la_stats->cat_clumps[*class], sizeof(double), compare);

		if(la_stats->cat_clumps[*class]%2!=0)
			median = list[la_stats->cat_clumps[*class]/2];
		else
			median = (list[la_stats->cat_clumps[*class]/2]+list[la_stats->cat_clumps[*class]/2-1])/2;
		free(list);
		return median;
	}
}

/* radius of gyration range */
double gyrate_ra(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double min=DBL_MAX;
	double max=0;

	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].gyrate>max)
				max = cl_stats[i].gyrate;
			if(cl_stats[i].gyrate<min)
				min = cl_stats[i].gyrate;
		}
		return max-min;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category!=*class)
				continue; // not this class
			if(cl_stats[i].gyrate>max)
				max = cl_stats[i].gyrate;
			if(cl_stats[i].gyrate<min)
				min = cl_stats[i].gyrate;
		}
		return max-min;
	}
}

/* radius of gyration standard deviation */
double gyrate_sd(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;

	if(!class){ // landscape level
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].gyrate;
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			dev = cl_stats[i].gyrate - mean;
			sd += dev*dev;
		}
		return sqrt(sd / la_stats->total_clumps);
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
		return 0.0; // class does not exist here
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].gyrate;
		mean = sum / la_stats->cat_clumps[*class];

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class){
				dev = cl_stats[i].gyrate - mean;
				sd += dev*dev;
			}
		}
		return sqrt(sd / la_stats->cat_clumps[*class]);
	}
}

/* radius of gyration coefficient of variation */
double gyrate_cv(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;

	if(!class){ // landscape level
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].gyrate;
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			dev = cl_stats[i].gyrate - mean;
			sd += dev*dev;
		}
		sd = sqrt(sd / la_stats->total_clumps);
		return sd / mean * 100;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].gyrate;
		mean = sum / la_stats->cat_clumps[*class];
		if(mean==0)
			return 0.0;

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class){
				dev = cl_stats[i].gyrate - mean;
				sd += dev*dev;
			}
		}
		sd = sqrt(sd / la_stats->cat_clumps[*class]);
		return sd / mean * 100;	
	}
}

/*============================Shape indices===================================*/

/* perimeter-area fractal dimension */
double pafrac(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double slope,nom,denom;
	double ln_p, ln_a;
	double sum_p=0,sum_a=0,sum_p2=0,sum_p_a=0;

	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			ln_p = log(cl_stats[i].perimeter);
			ln_a = log(cl_stats[i].area);
			sum_p += ln_p;
			sum_a += ln_a;
			sum_p2 += ln_p * ln_p;
			sum_p_a += ln_p * ln_a;
		}
		nom = (la_stats->total_clumps*sum_p_a)-(sum_p*sum_a);
		denom = (la_stats->total_clumps*sum_p2)-(sum_p*sum_p);
		if(!nom || !denom){
			return 1;
		}
		else{
			slope = nom / denom;
			return 2/slope;
		}
	}
	else{ // single class level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category!=*class)
				continue; // not this class
			ln_p = log(cl_stats[i].perimeter);
			ln_a = log(cl_stats[i].area);
			sum_p += ln_p;
			sum_a += ln_a;
			sum_p2 += ln_p * ln_p;
			sum_p_a += ln_p * ln_a;
		}
		nom = (la_stats->cat_clumps[*class]*sum_p_a)-(sum_p*sum_a);
		denom = (la_stats->cat_clumps[*class]*sum_p2)-(sum_p*sum_p);
		if(!nom || !denom){
			return 1;
		}
		else{
			slope = nom / denom;
			return 2/slope;
		}
	}
}

/* perimeter-area ratio */
/* perimeter-area ratio mean */
double para_mn(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double* para;
	
	if(!class){ // landscape level
		para=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			para[i-1] = cl_stats[i].perimeter / cl_stats[i].area;
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += para[i];
		free(para);
		return sum / la_stats->total_clumps;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		para=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				para[current++] = cl_stats[i].perimeter / cl_stats[i].area;
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += para[i];
		free(para);
		return sum / la_stats->cat_clumps[*class];
	}
}

/* perimeter-area ratio area-weighted mean */
double para_am(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i; 
	double sum=0;
	double* para;
	
	if(!class){ // landscape level
		para=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			para[i-1] = cl_stats[i].perimeter / cl_stats[i].area;
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += para[i] * (cl_stats[i+1].area / la_stats->notnull_area);
		free(para);
		return sum;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		para=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				para[current++] = cl_stats[i].perimeter / cl_stats[i].area;
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += para[i] * (cl_stats[i+1].area / la_stats->cat_areas[*class]);
		free(para);
		return sum;
	}
}

/* perimeter-area ratio median */
double para_md(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double median;
	double* para;
	
	if(!class){ // landscape level
		para=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			para[i-1] = cl_stats[i].perimeter / cl_stats[i].area;
		qsort(para, la_stats->total_clumps, sizeof(double), compare);

		if(la_stats->total_clumps%2)
			median = para[la_stats->total_clumps/2];
		else
			median = (para[la_stats->total_clumps/2]+para[la_stats->total_clumps/2-1])/2;
		free(para);
		return median;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		para=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				para[current++] = cl_stats[i].perimeter / cl_stats[i].area;

		if(la_stats->cat_clumps[*class]%2)
			median = para[la_stats->cat_clumps[*class]/2];
		else
			median = (para[la_stats->cat_clumps[*class]/2]+para[la_stats->cat_clumps[*class]/2-1])/2;
		free(para);
		return median;
	}
}

/* perimeter-area ratio range */
double para_ra(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double min=DBL_MAX;
	double max=0;
	double* para;
	
	if(!class){ // landscape level
		para=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			para[i-1] = cl_stats[i].perimeter / cl_stats[i].area;
		
		for(i=0; i<la_stats->total_clumps; ++i){
			if(para[i]>max)
				max = para[i];
			if(para[i]<min)
				min = para[i];
		}
		free(para);
		return max-min;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		para=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				para[current++] = cl_stats[i].perimeter / cl_stats[i].area;

		for(i=0; i<la_stats->cat_clumps[*class]; ++i){
			if(para[i]>max)
				max = para[i];
			if(para[i]<min)
				min = para[i];
		}
		free(para);
		return max-min;		
	}
}

/* perimeter-area ratio standard deviation */
double para_sd(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;
	double* para;
	
	if(!class){ // landscape level
		para=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			para[i-1] = cl_stats[i].perimeter / cl_stats[i].area;
		
		/* mean */
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += para[i];
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=0; i<la_stats->total_clumps; ++i){
			dev = para[i] - mean;
			sd += dev*dev;
		}
		free(para);
		return sqrt(sd / la_stats->total_clumps);
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
		return 0.0; // class does not exist here
		long int current=0;
		para=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				para[current++] = cl_stats[i].perimeter / cl_stats[i].area;

		/* mean */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += para[i];
		mean = sum / la_stats->cat_clumps[*class];

		/* sd */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i){
			dev = para[i] - mean;
			sd += dev*dev;
		}
		free(para);
		return sqrt(sd / la_stats->cat_clumps[*class]);
	}
}

/* perimeter-area ratio coefficient of variation */
double para_cv(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;
	double* para;
	
	if(!class){ // landscape level
		para=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			para[i-1] = cl_stats[i].perimeter / cl_stats[i].area;
		
		/* mean */
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += para[i];
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=0; i<la_stats->total_clumps; ++i){
			dev = para[i] - mean;
			sd += dev*dev;
		}
		sd = sqrt(sd / la_stats->total_clumps);
		free(para);
		return sd / mean * 100;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
		return 0.0; // class does not exist here
		long int current=0;
		para=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				para[current++] = cl_stats[i].perimeter / cl_stats[i].area;

		/* mean */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += para[i];
		mean = sum / la_stats->cat_clumps[*class];

		/* sd */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i){
			dev = para[i] - mean;
			sd += dev*dev;
		}
		sd = sqrt(sd / la_stats->cat_clumps[*class]);
		free(para);
		return sd / mean * 100;		
	}
}

/* shape index */
/* shape index mean */
double shape_mn(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double* shape;
	
	if(!class){ // landscape level
		shape=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			shape[i-1] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += shape[i];
		free(shape);
		return sum / la_stats->total_clumps;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		shape=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				shape[current++] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += shape[i];
		free(shape);
		return sum / la_stats->cat_clumps[*class];
	}
}

/* shape index area-weighted mean */
double shape_am(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double* shape;
	
	if(!class){ // landscape level
		shape=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			shape[i-1] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += shape[i] * (cl_stats[i+1].area / la_stats->notnull_area);
		free(shape);
		return sum;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		shape=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				shape[current++] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += shape[i] * (cl_stats[i+1].area / la_stats->cat_areas[*class]);
		free(shape);
		return sum;
	}
}

/* shape index median */
double shape_md(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double median;
	double* shape;
	
	if(!class){ // landscape level
		shape=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			shape[i-1] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);
		qsort(shape, la_stats->total_clumps, sizeof(double), compare);

		if(la_stats->total_clumps%2)
			median = shape[la_stats->total_clumps/2];
		else
			median = (shape[la_stats->total_clumps/2]+shape[la_stats->total_clumps/2-1])/2;
		free(shape);
		return median;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		shape=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				shape[current++] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);

		if(la_stats->cat_clumps[*class]%2)
			median = shape[la_stats->cat_clumps[*class]/2];
		else
			median = (shape[la_stats->cat_clumps[*class]/2]+shape[la_stats->cat_clumps[*class]/2-1])/2;
		free(shape);
		return median;
	}
}

/* shape index range */
double shape_ra(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double min=DBL_MAX;
	double max=0;
	double* shape;
	
	if(!class){ // landscape level
		shape=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			shape[i-1] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);
		
		for(i=0; i<la_stats->total_clumps; ++i){
			if(shape[i]>max)
				max = shape[i];
			if(shape[i]<min)
				min = shape[i];
		}
		free(shape);
		return max-min;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		shape=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				shape[current++] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);

		for(i=0; i<la_stats->cat_clumps[*class]; ++i){
			if(shape[i]>max)
				max = shape[i];
			if(shape[i]<min)
				min = shape[i];
		}
		free(shape);
		return max-min;		
	}
}

/* shape index standard deviation */
double shape_sd(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;
	double* shape;
	
	if(!class){ // landscape level
		shape=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			shape[i-1] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);
		
		/* mean */
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += shape[i];
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=0; i<la_stats->total_clumps; ++i){
			dev = shape[i] - mean;
			sd += dev*dev;
		}
		free(shape);
		return sqrt(sd / la_stats->total_clumps);
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
		return 0.0; // class does not exist here
		long int current=0;
		shape=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				shape[current++] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);

		/* mean */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += shape[i];
		mean = sum / la_stats->cat_clumps[*class];

		/* sd */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i){
			dev = shape[i] - mean;
			sd += dev*dev;
		}
		free(shape);
		return sqrt(sd / la_stats->cat_clumps[*class]);
	}
}

/* shape index coefficient of variation */
double shape_cv(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;
	double* shape;
	
	if(!class){ // landscape level
		shape=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			shape[i-1] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);
		
		/* mean */
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += shape[i];
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=0; i<la_stats->total_clumps; ++i){
			dev = shape[i] - mean;
			sd += dev*dev;
		}
		sd = sqrt(sd / la_stats->total_clumps);
		free(shape);
		return sd / mean * 100;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
		return 0.0; // class does not exist here
		long int current=0;
		shape=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				shape[current++] = 0.25 * cl_stats[i].perimeter / sqrt(cl_stats[i].area);

		/* mean */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += shape[i];
		mean = sum / la_stats->cat_clumps[*class];

		/* sd */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i){
			dev = shape[i] - mean;
			sd += dev*dev;
		}
		sd = sqrt(sd / la_stats->cat_clumps[*class]);
		free(shape);
		return sd / mean * 100;		
	}
}

/* fractal dimension index */
/* fractal dimension index mean */
double frac_mn(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double* frac;
	
	if(!class){ // landscape level
		frac=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			frac[i-1] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += frac[i];
		free(frac);
		return sum / la_stats->total_clumps;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		frac=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				frac[current++] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += frac[i];
		free(frac);
		return sum / la_stats->cat_clumps[*class];
	}
}

/* fractal dimension index area-weighted mean */
double frac_am(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double* frac;
	
	if(!class){ // landscape level
		frac=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			frac[i-1] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += frac[i] * (cl_stats[i+1].area / la_stats->notnull_area);
		free(frac);
		return sum;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		frac=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				frac[current++] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += frac[i] * (cl_stats[i+1].area / la_stats->cat_areas[*class]);
		free(frac);
		return sum;
	}
}

/* fractal dimension index median */
double frac_md(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double median;
	double* frac;
	
	if(!class){ // landscape level
		frac=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			frac[i-1] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);
		qsort(frac, la_stats->total_clumps, sizeof(double), compare);

		if(la_stats->total_clumps%2)
			median = frac[la_stats->total_clumps/2];
		else
			median = (frac[la_stats->total_clumps/2]+frac[la_stats->total_clumps/2-1])/2;
		free(frac);
		return median;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		frac=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				frac[current++] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);

		if(la_stats->cat_clumps[*class]%2)
			median = frac[la_stats->cat_clumps[*class]/2];
		else
			median = (frac[la_stats->cat_clumps[*class]/2]+frac[la_stats->cat_clumps[*class]/2-1])/2;
		free(frac);
		return median;
	}
}

/* fractal dimension index range */
double frac_ra(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double min=DBL_MAX;
	double max=0;
	double* frac;
	
	if(!class){ // landscape level
		frac=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			frac[i-1] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);
		
		for(i=0; i<la_stats->total_clumps; ++i){
			if(frac[i]>max)
				max = frac[i];
			if(frac[i]<min)
				min = frac[i];
		}
		free(frac);
		return max-min;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		long int current=0;
		frac=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				frac[current++] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);

		for(i=0; i<la_stats->cat_clumps[*class]; ++i){
			if(frac[i]>max)
				max = frac[i];
			if(frac[i]<min)
				min = frac[i];
		}
		free(frac);
		return max-min;		
	}
}

/* fractal dimension index standard deviation */
double frac_sd(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;
	double* frac;
	
	if(!class){ // landscape level
		frac=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			frac[i-1] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);
		
		/* mean */
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += frac[i];
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=0; i<la_stats->total_clumps; ++i){
			dev = frac[i] - mean;
			sd += dev*dev;
		}
		free(frac);
		return sqrt(sd / la_stats->total_clumps);
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
		return 0.0; // class does not exist here
		long int current=0;
		frac=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				frac[current++] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);

		/* mean */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += frac[i];
		mean = sum / la_stats->cat_clumps[*class];

		/* sd */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i){
			dev = frac[i] - mean;
			sd += dev*dev;
		}
		free(frac);
		return sqrt(sd / la_stats->cat_clumps[*class]);
	}
}

/* fractal dimension index coefficient of variation */
double frac_cv(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;
	double* frac;
	
	if(!class){ // landscape level
		frac=malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			frac[i-1] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);
		
		/* mean */
		for(i=0; i<la_stats->total_clumps; ++i)
			sum += frac[i];
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=0; i<la_stats->total_clumps; ++i){
			dev = frac[i] - mean;
			sd += dev*dev;
		}
		sd = sqrt(sd / la_stats->total_clumps);
		free(frac);
		return sd / mean * 100;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
		return 0.0; // class does not exist here
		long int current=0;
		frac=calloc(la_stats->cat_clumps[*class],sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				frac[current++] = 2 * log(0.25*cl_stats[i].perimeter) / log(cl_stats[i].area);

		/* mean */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i)
			sum += frac[i];
		mean = sum / la_stats->cat_clumps[*class];

		/* sd */
		for(i=0; i<la_stats->cat_clumps[*class]; ++i){
			dev = frac[i] - mean;
			sd += dev*dev;
		}
		sd = sqrt(sd / la_stats->cat_clumps[*class]);
		free(frac);
		return sd / mean * 100;		
	}
}

/* contiguity index */
/* contiguity index mean */
double contig_mn(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].contig;
		return sum / la_stats->total_clumps;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].contig;
		return (sum / la_stats->cat_clumps[*class]) / 10000;
	}
}

/* contiguity index area-weighted mean */
double contig_am(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].contig * (cl_stats[i].area / la_stats->notnull_area);
		return sum;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += (cl_stats[i].contig * (cl_stats[i].area / la_stats->cat_areas[*class])) / 10000;
		return sum;
	}
}

/* contiguity index median */
double contig_md(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double median;
	double* list;

	if(!class){ // landscape level
		/* create temporary table and quicksort it */
		list = malloc(la_stats->total_clumps*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			list[i-1] = cl_stats[i].contig;
		qsort(list, la_stats->total_clumps, sizeof(double), compare);

		if(la_stats->total_clumps%2)
			median = list[la_stats->total_clumps/2];
		else
			median = (list[la_stats->total_clumps/2]+list[la_stats->total_clumps/2-1])/2;
		free(list);
		return median;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		/* create temporary table and quicksort it */
		long int current=0;
		list = malloc(la_stats->cat_clumps[*class]*sizeof(double));
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				list[current++] = cl_stats[i].contig;
		qsort(list, la_stats->cat_clumps[*class], sizeof(double), compare);

		if(la_stats->cat_clumps[*class]%2!=0)
			median = list[la_stats->cat_clumps[*class]/2];
		else
			median = (list[la_stats->cat_clumps[*class]/2]+list[la_stats->cat_clumps[*class]/2-1])/2;
		free(list);
		return median;
	}
}

/* contiguity index range */
double contig_ra(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double min=DBL_MAX;
	double max=0;

	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].contig>max)
				max = cl_stats[i].contig;
			if(cl_stats[i].contig<min)
				min = cl_stats[i].contig;
		}
		return max-min;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category!=*class)
				continue; // not this class
			if(cl_stats[i].contig>max)
				max = cl_stats[i].contig;
			if(cl_stats[i].contig<min)
				min = cl_stats[i].contig;
		}
		return max-min;
	}
}

/* contiguity index standard deviation */
double contig_sd(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;

	if(!class){ // landscape level
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].contig;
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			dev = cl_stats[i].contig - mean;
			sd += dev*dev;
		}
		return sqrt(sd / la_stats->total_clumps);
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
		return 0.0; // class does not exist here
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].contig;
		mean = sum / la_stats->cat_clumps[*class];

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class){
				dev = cl_stats[i].contig - mean;
				sd += dev*dev;
			}
		}
		return sqrt(sd / la_stats->cat_clumps[*class]);
	}
}

/* contiguity index coefficient of variation */
double contig_cv(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	double mean,dev,sd=0;

	if(!class){ // landscape level
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			sum += cl_stats[i].contig;
		mean = sum / la_stats->total_clumps;

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			dev = cl_stats[i].contig - mean;
			sd += dev*dev;
		}
		sd = sqrt(sd / la_stats->total_clumps);
		return sd / mean * 100;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		/* mean */
		for(i=1; i<la_stats->total_clumps+1; ++i)
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].contig;
		mean = sum / la_stats->cat_clumps[*class];
		if(mean==0)
			return 0.0;

		/* sd */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class){
				dev = cl_stats[i].contig - mean;
				sd += dev*dev;
			}
		}
		sd = sqrt(sd / la_stats->cat_clumps[*class]);
		return sd / mean * 100;	
	}
}

/*============================Aggregation indices=============================*/

/* contagion index (%) */
double contag(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	/* function uses "double-count" adjacency */
	/* if there is only one category in the landscape, contag can't be calculated. 
	In such a case, FRAGSTATS returns N/A.
	However, for further similarity calculation, we need to have a number in the
	output feature vector. In teory, 100 is a correct contagion value, however, 
	only if the single patch covers the entire landscape. If there exist nulls 
	in the landscape, it will result in a weaker inverse correlation between 
	contagion and edge density (which are naturally strongly correlated).
	In other words: the function calculates incorrect value if there is only 
	one category in the landscape, and there are null values in the landscape. */
	if(!class){ // landscape level
		if(la_stats->ncats < 2){
			return 100;
		}
		int i,k;
		long int sum_g,adj,cat_i;
		double tmp,sum=0;
		int length = map_cats->cat[map_cats->num-1] + 1; // highest category number + 1
		double prop_i; // category proportion of landscape

		/* for each category */
		for(i=0; i<map_cats->num; ++i){
			cat_i = map_cats->cat[i];
			prop_i = la_stats->cat_areas[cat_i] / la_stats->notnull_area;

			/* sum up adjacencies of category i (double count) */
			sum_g=0;
			for(k=0; k<map_cats->num; ++k)
				sum_g += la_stats->adjacency[cat_i*length+map_cats->cat[k]];

			if(sum_g==0) // category i has no notnull neighbors or doesn't appear in the landscape
				continue;
			/* sum over all the adjacencies */
			for(k=0; k<map_cats->num; ++k){
				if((adj=la_stats->adjacency[cat_i*length+map_cats->cat[k]]) == 0)
					continue; // no adjacencies for the i-k combination
				tmp = prop_i * ((double)adj / sum_g);
				sum += tmp * log(tmp);
			}
		} // end for i
		return ( 1 + sum/(2*log(la_stats->ncats)) ) * 100;
	}
	else{ // single class level
		G_fatal_error("Index can be calculated only for landscape level.");
	}
	return 0;
}

/* interspersion & juxtaposition index (%) */
double iji(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	/* if there is less than 3 categories in the landscape, IJI is undefined. 
	In such a case, FRAGSTATS returns N/A.
	However, for further similarity calculation, we need to have a number in the
	output feature vector. 100 is set (meaning all patch types are equally 
	adjacent to each other). */
	if(la_stats->ncats < 3){
		return 100;
	}

	int i,k;
	int length = map_cats->cat[map_cats->num-1] + 1; // highest category number + 1
	double prop,edge,total_edge=0,sum=0;

	if(!class){ // landscape level
		/* sum of notnull edges */
		for(i=0; i<map_cats->num; ++i)
			for(k=i+1; k<map_cats->num; ++k)
				total_edge += la_stats->cat_edge_matrix[map_cats->cat[i]*length+map_cats->cat[k]];

		/* for each category */
		for(i=0; i<map_cats->num; ++i){
			for(k=i+1; k<map_cats->num; ++k){
				if((edge=la_stats->cat_edge_matrix[map_cats->cat[i]*length+map_cats->cat[k]]) == 0)
					continue; // classes i and k are not adjacent
				prop = edge / total_edge;
				sum += prop * log(prop);
			}
		}
		return -sum / (log(0.5*la_stats->ncats*(la_stats->ncats-1))) * 100;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		/* sum of notnull edges */
		for(k=0; k<map_cats->num; ++k)
			if(*class!=map_cats->cat[k])
				total_edge += la_stats->cat_edge_matrix[(*class)*length+map_cats->cat[k]];
		
		for(k=0; k<map_cats->num; ++k){
			if((edge=la_stats->cat_edge_matrix[(*class)*length+map_cats->cat[k]]) == 0)
				continue; // classes i and k are not adjacent
			prop = edge / total_edge;
			sum += prop * log(prop);
		}
		return -sum / (log(la_stats->ncats-1)) * 100;
	}
}

/* percentage of like adjacencies (%) */
double pladj(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	/* FRAGSTATS returns N/A if the landscape consists of a single not-null	cell.
	However, mathematically, index is defined for this case and it is 0. */

	int i,cat;
	long int sum=0;
	int length = map_cats->cat[map_cats->num-1] + 1; // highest category number + 1
	
	if(!class){ // landscape level
		/* sum up like adjacencies */
		for(i=0; i<map_cats->num; ++i){
			cat = map_cats->cat[i];
			sum += la_stats->adjacency[cat*length+cat];
		}
		return ((double)sum / la_stats->total_adj) * 100;
	}
	else{ // single class level
		if(la_stats->cat_clumps[*class]==0)
			return 0.0; // class does not exist here
		return ((double)la_stats->adjacency[(*class)*length+(*class)] / la_stats->cat_adj[*class]) * 100;
	}
}

/* aggregation index (%) */
double ai(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	/* Index is undefined if each existing category consists of a single cell,
	in which case FRAGSTATS returns N/A.
	However, for further similarity calculation, we need to have a number in the
	output feature vector. 0 is set (meaning max disaggregation of categories). */

	int i,cat,sth_done=0;
	int length = map_cats->cat[map_cats->num-1] + 1; // highest category number + 1
	double res_area = res*res; // map areal resolution
	long int Ai,maxi,ni,mi;
	double sum=0;
	double prop; // category proportion of landscape

	if(!class){ // landscape level
		/* sum over all categories */
		for(i=0; i<map_cats->num; ++i){
			cat = map_cats->cat[i];
			prop = la_stats->cat_areas[cat] / la_stats->notnull_area;
			if((Ai = (long)(la_stats->cat_areas[cat]/res_area)) < 2)
				continue; // too few cells
			ni=(long)sqrt(Ai);
			mi=Ai - ni*ni;
			if(mi==0)
				maxi=2*ni*(ni-1);
			else if(mi>ni)
				maxi=2*ni*(ni-1) + 2*mi-2;
			else
				maxi=2*ni*(ni-1) + 2*mi-1;
			/* single-count */
			sum+= (double)(la_stats->adjacency[cat*length+cat]/2)/maxi * prop;
			sth_done=1;
		}

		if(sth_done)
			return sum * 100;
		else // there is no category with cell_count>1
			return 0;
	}
	else{ // single class level
		if((Ai = (long)(la_stats->cat_areas[*class]/res_area)) < 2)
			return 0; // too few cells
		ni=(long)sqrt(Ai);
		mi=Ai - ni*ni;
		if(mi==0)
			maxi=2*ni*(ni-1);
		else if(mi>ni)
			maxi=2*ni*(ni-1) + 2*mi-2;
		else
			maxi=2*ni*(ni-1) + 2*mi-1;	
		/* single-count */
		return (double)(la_stats->adjacency[(*class)*length+(*class)]/2)/maxi * 100;
	}
}

/* landscape shape index */
double lsi(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	if(!class){ // landscape level
		/* add boundary to the total edge */
		double tedge = la_stats->total_edge + 2*(nrows+ncols)*res;
		return 0.25*tedge / sqrt(la_stats->notnull_area);
	}
	else{ // single class level
		unsigned long int i;
		double tedge=0;
		/* sum up perimeters of patches to get sum of all edges for the class */
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class)
				tedge += cl_stats[i].perimeter;
		}
		return 0.25*tedge / sqrt(la_stats->notnull_area);
	}
}

/* patch cohesion index */
double cohesion(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	unsigned int perim_cells; // perimeter in cells
	unsigned int total_cells = (unsigned int)(la_stats->notnull_area / (res*res));
	double sum_p=0, sum_pa=0;
	double res_area = res*res; // map areal resolution

	if(!class){ // landscape level
		if((la_stats->notnull_area/res_area) < 2)
			return 0; // too few cells
		for(i=1; i<la_stats->total_clumps+1; ++i){
			perim_cells = (unsigned int)(cl_stats[i].perimeter / res);
			sum_p += perim_cells;
			sum_pa += perim_cells * sqrt(cl_stats[i].cells);
		}
		return (1-sum_p/sum_pa)/(1-1/sqrt(total_cells))*100;
	}
	else{ // single class level
		if((la_stats->cat_areas[*class]/res_area) < 2)
			return 0; // too few cells
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class){
				perim_cells = (unsigned int)(cl_stats[i].perimeter / res);
				sum_p += perim_cells;
				sum_pa += perim_cells * sqrt(cl_stats[i].cells);
			}
		}
		return (1-sum_p/sum_pa)/(1-1/sqrt(total_cells))*100;
	}
}

/* patch density (amount/100ha) */
double pd(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	if(!class){ // landscape level
		return (double)la_stats->total_clumps/la_stats->notnull_area * 1000000;
	}
	else{ // single class level
		return (double)la_stats->cat_clumps[*class]/la_stats->notnull_area * 1000000;
	}
}

/* splitting index */
double split(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	
	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			sum += cl_stats[i].area * cl_stats[i].area;
		}
		return la_stats->notnull_area*la_stats->notnull_area / sum;
	}
	else{ // single class level
		if(!la_stats->cat_clumps[*class])
			return 0; // class does not exist
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].area * cl_stats[i].area;
		}
		return la_stats->notnull_area*la_stats->notnull_area / sum;
	}
}

/* landscape division index (proportion) */
double division(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double prop;
	double sum=0;
	
	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			prop = cl_stats[i].area/la_stats->notnull_area;
			sum += prop*prop;
		}
		return 1- sum;
	}
	else{ // single class level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class){
				prop = cl_stats[i].area/la_stats->notnull_area;
				sum += prop*prop;
			}
		}
		return 1- sum;
	}
}

/* effective mesh size (ha) */
double mesh(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	unsigned long int i;
	double sum=0;
	
	if(!class){ // landscape level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			sum += cl_stats[i].area * cl_stats[i].area;
		}
		return (sum / la_stats->notnull_area) / 10000 ;
	}
	else{ // single class level
		for(i=1; i<la_stats->total_clumps+1; ++i){
			if(cl_stats[i].category==*class)
				sum += cl_stats[i].area * cl_stats[i].area;
		}
		return (sum / la_stats->notnull_area) / 10000 ;	
	}
}


/*=============================Diversity indices==============================*/

/* patch richness density (amount/100ha) */
double prd(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	if(!class){ // landscape level
		return (double)la_stats->ncats / la_stats->notnull_area * 1000000;
	}
	else{ // single class level
		G_fatal_error("Index can be calculated only for landscape level.");
	}
	return 0;
}

/* relative patch richness (%) */
double rpr(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	if(!class){ // landscape level
		return (double)la_stats->ncats / map_cats->num * 100;
	}
	else{ // single class level
		G_fatal_error("Index can be calculated only for landscape level.");
	}
	return 0;
}

/* Shannon's diversity index */
double shdi(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	long int i;
	double prop,sum=0,area;
	
	if(!class){ // landscape level
		for(i=0; i<map_cats->num; ++i){
			if((area=la_stats->cat_areas[map_cats->cat[i]]) == 0)
				continue; // class does not appear in the landscape
			prop = area / la_stats->notnull_area;
			sum += prop * log(prop);
		}
		return -sum;
	}
	else{ // single class level
		G_fatal_error("Index can be calculated only for landscape level.");
	}
	return 0;
}

/* Simpson's diversity index */
double sidi(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	long int i;
	double prop,sum=0;
	
	if(!class){ // landscape level
		for(i=0; i<map_cats->num; ++i){
			prop = la_stats->cat_areas[map_cats->cat[i]] / la_stats->notnull_area;
			sum += prop * prop;
		}
		return 1 - sum;
	}
	else{ // single class level
		G_fatal_error("Index can be calculated only for landscape level.");
	}
	return 0;
}

/* modified Simpson's diversity index */
double msidi(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	long int i;
	double prop,sum=0;
	
	if(!class){ // landscape level
		for(i=0; i<map_cats->num; ++i){
			prop = la_stats->cat_areas[map_cats->cat[i]] / la_stats->notnull_area;
			sum += prop * prop;
		}
		return -log(sum);
	}
	else{ // single class level
		G_fatal_error("Index can be calculated only for landscape level.");
	}
	return 0;
}

/* Shannon's evenness index */
double shei(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	long int i;
	double prop,sum=0,area;
	
	if(!class){ // landscape level
		if(la_stats->ncats<2){
			return 0;
		}
	
		for(i=0; i<map_cats->num; ++i){
			if((area=la_stats->cat_areas[map_cats->cat[i]]) == 0)
				continue; // class does not appear in the landscape
			prop = area / la_stats->notnull_area;
			sum += prop * log(prop);
		}
		return (-sum) / log(la_stats->ncats);
	}
	else{ // single class level
		G_fatal_error("Index can be calculated only for landscape level.");
	}
	return 0;
}

/*  Simpson's evenness index */
double siei(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	long int i;
	double prop,sum=0;

	if(!class){ // landscape level
		if(la_stats->ncats<2){
			return 0;
		}

		for(i=0; i<map_cats->num; ++i){
			prop = la_stats->cat_areas[map_cats->cat[i]] / la_stats->notnull_area;
			sum += prop * prop;
		}
		return (1.-sum) / (1.-1./la_stats->ncats);
	}
	else{ // single class level
		G_fatal_error("Index can be calculated only for landscape level.");
	}
	return 0;
}

/* modified Simpson's evenness index */
double msiei(EZGDAL_LAYER* l, LI_CLUMP* cl_stats, LI_LANDSCAPE* la_stats,
	long int nrows, long int ncols, double res, MAP_CATS* map_cats, int* class)
{
	long int i;
	double prop,sum=0;

	if(!class){ // landscape level
		if(la_stats->ncats<2){
			return 0;
		}

		for(i=0; i<map_cats->num; ++i){
			prop = la_stats->cat_areas[map_cats->cat[i]] / la_stats->notnull_area;
			sum += prop * prop;
		}
		return (-log(sum)) / log(la_stats->ncats);
	}
	else{ // single class level
		G_fatal_error("Index can be calculated only for landscape level.");
	}
	return 0;
}

