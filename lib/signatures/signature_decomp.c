/****************************************************************************
 *
 * MODULE:	simple decomposition signature
 * AUTHOR(S):	Jacek Niesterowicz, Jaroslaw Jasiewicz
 * PURPOSE:	calculating decomposition signature:
 *		simple version
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <stdlib.h>
#include <ezgdal.h>
#include <stdarg.h>

#include "signature_geopat_compatibility.h"
#include "signature_decomposition.h"

/* SIMPLE 2-level DECOMPOSITION
 * decomposition at 2 levels: 4x4 and 16x16;
 * Parameters can be changed at the beginning of the main function
 * If parameters are changed, change numbers in decomposition_len function.
 * NOTE: decomposition requires square-shaped input pattern */

int decompose_area(EZGDAL_LAYER* layer, EZGDAL_FRAME* f, DC_PARAMS* p, int row, int col, double *signature)
{
	int category,fine_index,index;
	int nulls=0;
	int* th_coarse;
	int* th_fine;
	int coarse_class, fine_class;
	int decomp_level=p->dc_level;
	int r,c,i,j;
	int er, ec;
	int num_of_fine_quads=decomp_level*decomp_level;
	int size_of_fine_quad=p->dc_base_size;
	int area_of_fine_quad=size_of_fine_quad*size_of_fine_quad;
	int num_of_classes=(int)log2(p->dc_region_size*p->dc_region_size);
	int ncats = layer->stats->map_max_val+1;
	
	th_coarse=calloc(ncats,sizeof(int));
	th_fine=calloc(ncats*num_of_fine_quads,sizeof(int)); /* region is divided into 4*4 =16 subregions */

	for(r=0;r<p->dc_region_size;++r)
		for(c=0;c<p->dc_region_size;++c) {
			er=r+row;
			ec=c+col;
			if(ezgdal_is_null(layer, f->buffer[er][ec])){
				nulls+=1;
				continue;
			}
			category=ezgdal_get_value_index(layer, f->buffer[er][ec]);
			th_coarse[category]++;

			fine_index=(r/size_of_fine_quad)*decomp_level+(c/size_of_fine_quad); /* integer division to determine sub-region */
			fine_index=fine_index*ncats+category;
			th_fine[fine_index]++;
		}

	for(i=0;i<ncats;++i){
		coarse_class=th_coarse[i]>area_of_fine_quad?(int)log2(th_coarse[i]-1):(int)log2(area_of_fine_quad);
		index=i*num_of_classes+coarse_class;
		signature[index]+=th_coarse[i];

		for(j=0;j<num_of_fine_quads;++j) { /* fine quads starts from 0 to 15 */
			fine_class=th_fine[j*ncats+i]>1?(int)log2(th_fine[j*ncats+i]-1):0;
			index=i*num_of_classes+fine_class;
			signature[index]+=th_fine[j*ncats+i];
		}
	}

	free(th_coarse);
	free(th_fine);
	
	return nulls;
}

/* main function */
int decomposition(EZGDAL_FRAME **frames, int num_of_frames, double *signature, int signature_len, ...)
{
	int nrows, ncols;
	int r,c;
	int step,count=0;
	int nulls=0,total_cells;
	EZGDAL_FRAME *f = frames[0];
	EZGDAL_LAYER *l = f->owner.stripe->layer;
	DC_PARAMS *p = malloc(sizeof(DC_PARAMS));

	/* set parameters - can be modified */
	p->dc_base_size=4; //size of smallest quad (one dimension) area is 4*4
	p->dc_level=4; //decomposition level, fixed
	p->dc_region_size=p->dc_base_size*p->dc_level; // quad size at the lowest decomposition level

	if(f->rows!=f->cols)
		G_fatal_error("Decomposition: input pattern is not square (%d x %d)", f->rows, f->cols);

	if(p->dc_region_size>f->rows)
		G_fatal_error("Decomposition: pattern size must be at least %d", p->dc_region_size);

	/* lowering step enables overlaping */
	step=p->dc_region_size;
	nrows = f->rows-p->dc_region_size+1;
	ncols = f->cols-p->dc_region_size+1;

	for(r=0;r<nrows;r+=step)
		for(c=0;c<ncols;c+=step){
			nulls+=decompose_area(l,f,p,r,c,signature);
			count++;
		}

	total_cells=count*p->dc_region_size*p->dc_region_size;
	free(p);

	if(nulls*2 > total_cells)
		return 0; /* mark as null landscape */
	else
		return 1;
}

int decomposition_len(EZGDAL_LAYER **layers, int num_of_layers) {
 /* decomposition uses two dimentional histogram
  * one dimention are categories used in the map
  * second dimention is decomposition
  * first half of decomposition component is fine elements second is large element */
	int ncats = layers[0]->stats->map_max_val+1;
	/* hard-coded for 16x16 coarse size */
	return ((int)log2(16*16))*ncats;
}

