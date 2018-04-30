/****************************************************************************
 *
 * MODULE:	landscape indices signature
 * AUTHOR(S):	Jacek Niesterowicz
 * PURPOSE:	calculating a vector of landscape indices:
 *		full vector size
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include "signature_landind_macro.h"
#include "signature_landind_lips.h"
#include "../../lib/ezGDAL/ezgdal.h"
#include <stdarg.h>

int landind(EZGDAL_FRAME **frames, int num_of_frames, double *signature, int signature_len, ...)
{
	int i,j;
	long int row,col,r,c,next_r,next_c,nrows,ncols;
	int start,increment;
	long int index,target;
	long int first,last;
	int pos;
	int cat, target_cat;
	unsigned long int curr_clump=0;
	double** map; /* input map */
	unsigned long int* map_clump; /* clump map */
	long int* queue = NULL;
	LI_LANDSCAPE land_stats;
	LI_CLUMP* clump_stats=NULL;
	calculate_index* calculate;
	EZGDAL_FRAME *f = frames[0];
	EZGDAL_LAYER *l = f->owner.stripe->layer;
	H_PARAMS* p=malloc(sizeof(H_PARAMS));
	int length = li_init_landscape(l, &land_stats);
	MAP_CATS map_cats;

	/* get resolution */
	double* at = ezgdal_layer_get_at(l);
	double resolution = at[1]; // map resolution
	double res_area = resolution*resolution; // map areal resolution
	free(at);

	/* read map categories */
	map_cats.num = l->stats->map_max_val+1;
	map_cats.cat=malloc(map_cats.num*sizeof(int));
	j=0;
	for(i=0; i<(int)l->stats->hist_N; i++)
		if(l->stats->map_cat[i] >= 0)
			map_cats.cat[j++] = i + (int)(l->stats->min);

	/* SET HARDCODED PARAMETERS FOR INDICES */
	/* 0 - 4-neighborhood; 1 - 8-neighborhood */
	li_set_params_all(p,l,1);

	/* get size of the window */
	nrows = f->rows;
	ncols = f->cols;

	/* set easy-to-use input map buffer */
	map = malloc(nrows * sizeof(double*));
	for(row=0; row<nrows; ++row)
		map[row] = &(f->buffer[row][0]);

	/* allocate clump map */
	map_clump = calloc(nrows*ncols, sizeof(unsigned long int)); // calloc initializes to 0
	queue = malloc(4*nrows*ncols*sizeof(long int));

	/* prepare parameters for neighborhood-check loop */
	if(p->eight_flag){
		start=1;
		increment=1;
	}
	else{
		start=2;
		increment=2;
	}
/*printf("begin clump\n");*/
	/* Phase 1: clump & update statistics */
	for(row=0; row<nrows; ++row){
		for(col=0; col<ncols; ++col){
			land_stats.total_area+=res_area;
			if(map_clump[row*ncols+col] || ezgdal_is_null(l, map[row][col]))
				continue; // already clumped or null value or not in the circle

			curr_clump++;
			last=1;
			first=0;
			queue[0]=INDEX(row,col);
			clump_stats = li_add_clump(clump_stats, &land_stats, curr_clump, (int)map[row][col], row, col);
			cat = clump_stats[curr_clump].category;

			do {
				index=queue[first++];
				if(map_clump[index])
					continue;
				r=index/ncols;
				c=index%ncols;
				/* update clump and landscape area statistics */
				map_clump[index] = curr_clump;
				li_update_clump(clump_stats, &land_stats, curr_clump, (int)map[r][c], r, c, res_area, p->eight_flag);

				/* check neighbors */
				for(i=start;i<9;i+=increment) {
					next_r=SNR(i);
					next_c=SNC(i);

					if(NOT_IN_REGION(i)){
						if(IS_FOUR_CONN(i))
							clump_stats[curr_clump].perimeter+=resolution;
						continue;
					}
					target_cat = (int)map[next_r][next_c];
					
					if(ezgdal_is_null(l, (double)target_cat)) {
						if(IS_FOUR_CONN(i)){
							clump_stats[curr_clump].perimeter+=resolution;
							land_stats.total_edge+=resolution;
							land_stats.cat_edges[cat]+=resolution;
						}
						continue;
					}

					target=INDEX(next_r,next_c);

					/* collect adjacency and edge info */
					if(IS_FOUR_CONN(i)){
						land_stats.adjacency[cat*length+target_cat]++;
					}
 					if(cat==target_cat){
 						if(map_clump[target] == 0)
							queue[last++]=target;
					}
					else if(IS_FOUR_CONN(i)){
						clump_stats[curr_clump].perimeter+=resolution;
						land_stats.cat_edges[cat]+=resolution;
						land_stats.cat_edge_matrix[cat*length+target_cat]+=resolution;
						if(map_clump[target] == 0) // for total: don't count between-class edges twice
							land_stats.total_edge+=resolution;
					}
				}
			} while(first!=last);
		} // end col
	} // end row
/*printf("Calculate indices.\n");*/
	/* Phase 2: calculate post-clumping statistics */
	if(li_check_gyrate(p))
		li_gyrate(map_clump, clump_stats, land_stats.total_clumps, nrows, ncols, resolution);
	if(li_check_contig(p))
		li_contig(map_clump, clump_stats, land_stats.total_clumps, nrows, ncols, resolution);
	
	if(land_stats.notnull_area==0){ /* no data for the landscape */
		free(map);
		free(map_clump);
		free(queue);
		free(clump_stats);
		free(map_cats.cat);
		li_free_landscape(&land_stats);
		li_free_parameters(p);
		return 0;
	}

	pos=0; /* vector index */

	/* Phase 3: calculate landscape-level indices */
	if(p->li_land_names){
		for(i=0; i<li_indices_number(ALL); ++i){
			if(!p->li_land_flags[i])
				continue; // index marked as "don't calculate"

			/* set function and calculate the index */
			calculate = get_index_func(i);
			signature[pos++] = calculate(l, clump_stats, &land_stats, nrows, ncols, resolution, &map_cats, NULL);
		} // end for
	}

	/* Phase 4: calculate class-level indices */
	if(p->li_class_names){
		for(i=0; i<li_indices_number(ALL); ++i){
			if(!p->li_class_flags[i])
				continue; // index marked as "don't calculate"

			/* for each class */
			for(j=0; j<p->num_of_classes; ++j){
				/* set function and calculate the index */
				calculate = get_index_func(i);
				signature[pos++] = calculate(l, clump_stats, &land_stats, nrows, ncols, resolution, &map_cats, &p->classes[j]);
			}
		}
	}

	free(map);
	free(map_clump);
	free(queue);
	free(clump_stats);
	free(map_cats.cat);
	li_free_landscape(&land_stats);
	/* FREE HARDCODED PARAMETERS */
	li_free_parameters(p);

	if((land_stats.notnull_area/land_stats.total_area) < 0.5)
		return 0; /* mark as a null landscape */
	else
		return 1;
}

/* length of resulting vectors if ALL the indices are to be calculated */
int landind_len(EZGDAL_LAYER **layers, int num_of_layers)
{
	return li_indices_number(LANDSCAPE)+(li_indices_number(CLASS)*(layers[0]->stats->map_max_val+1));
}

