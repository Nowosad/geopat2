/****************************************************************************
 *
 * MODULE:	full decomposition signature
 * AUTHOR(S):	Jacek Niesterowicz, Jaroslaw Jasiewicz
 * PURPOSE:	calculating decomposition signature:
 *		full version
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <stdlib.h>
#include "../../lib/ezGDAL/ezgdal.h"
#include <math.h>
#include <stdarg.h>

#include "signature_geopat_compatibility.h"
#include "signature_decomposition.h"

/* FULL DECOMPOSITION
 * Additional parameter "dc_level" - number of decomposition levels.
 * size at the most coarse level = 2^(dc_level+1)
 * Default dc_level=0 - all possible levels for given pattern size
 * NOTE: decomposition requires square-shaped input pattern. */

int isPowerOfTwo(int x)
{
 while(((x%2)==0) && x>1) /* While x is even and > 1 */
   x/=2;
 return (x==1);
}

int get_size_division(int num, int quad_area, int num_of_size_divisions)
{
	int i;
	for(i=1;i<num_of_size_divisions;++i)
		if(num>(quad_area/(int)pow(2,i)))
			return num_of_size_divisions-i;
	return num_of_size_divisions-i;
}

int full_decompose_area(EZGDAL_LAYER* layer, EZGDAL_FRAME* f, DC_PARAMS* p, int row, int col, double *signature)
{
	int q,k,l,r,c,er,ec;
	int category,cat_index,quad_index,index;
	int num_of_quads_at_level;
	int total_quads_at_level;
	int quad_size_at_level;
	int quad_area_at_level;
//	int quad_weight;
	int total_quads=p->dc_num_of_quads*p->dc_num_of_quads;
	int ncats = layer->stats->map_max_val+1;
	int nulls=0;
	int** level_quads=malloc(p->dc_level*sizeof(int*));
	
	for(l=0;l<p->dc_level;++l) //allocate memory for level_QUADS and init to 0
		level_quads[l]=calloc((int)(total_quads/pow(4,l))*ncats,sizeof(int));

	/* scan map */
	for(r=0;r<p->dc_region_size;++r)
		for(c=0;c<p->dc_region_size;++c) {
			er=r+row;
			ec=c+col;
			category = (int)(f->buffer[er][ec]);
			if(ezgdal_is_null(layer, f->buffer[er][ec])){
				nulls++;
				continue;
			}
			cat_index=ezgdal_get_value_index(layer, category);
			/* update quad's histogram at each level */
			for(l=0;l<p->dc_level;++l){
				quad_size_at_level=p->dc_base_size*(int)pow(2,l);
				num_of_quads_at_level=p->dc_num_of_quads/(int)pow(2,l);
				quad_index=(r/quad_size_at_level)*num_of_quads_at_level+(c/quad_size_at_level);
				level_quads[l][quad_index*ncats+cat_index]++;
			}
	} //end for

	/* update actual full_decomposition histograms */
	for(l=0;l<p->dc_level;++l) {
		total_quads_at_level=total_quads/(int)pow(4,l);
		quad_size_at_level=p->dc_base_size*(int)pow(2,l);
		quad_area_at_level=quad_size_at_level*quad_size_at_level;
//		quad_weight=pow(4,l); // equals to the number of base quads

		for(q=0;q<total_quads_at_level;++q) {
			for(k=0;k<ncats;++k) {
				index=l*ncats*p->dc_num_of_size_divisions+k*p->dc_num_of_size_divisions+
					get_size_division(level_quads[l][q*ncats+k],quad_area_at_level,p->dc_num_of_size_divisions);
//				signature[index]+=quad_weight;
				signature[index]+=level_quads[l][q*ncats+k];
			}
		}
	}

	for(l=0;l<p->dc_level;++l)
		free(level_quads[l]);
	free(level_quads);

	return nulls;
}

/* main function */
int full_decomposition(EZGDAL_FRAME **frames, int num_of_frames, double *signature, int signature_len, ...)
{

	int nrows, ncols;
	int r,c;
	int step;
	int min_size;
	int nulls=0,total_cells;
	int count=0;
	int pattern_size;
	EZGDAL_FRAME *f = frames[0];
	EZGDAL_LAYER *l = f->owner.stripe->layer;
	DC_PARAMS *p = malloc(sizeof(DC_PARAMS));

	int dc_level;
	va_list arg_list;
	va_start(arg_list,signature_len);
	dc_level = va_arg(arg_list, int);
	va_end(arg_list);


	if(f->rows!=f->cols)
		G_fatal_error("Full decomposition: input pattern is not square (%d x %d)", f->rows, f->cols);

	
	/* set full_decomp parameters - can be modified */

	p->dc_num_of_size_divisions=3; //internal constant
	p->dc_base_size=4; //internal constant (must be power of 2)
	min_size=p->dc_base_size; //internal constant

	/* get and check pattern_size */

	pattern_size = f->rows;
	if(pattern_size<min_size)
		G_fatal_error("Full decomposition: pattern size must be at least %d.",min_size);

	/* custom or dynamic set of the most coarse level */

	if(dc_level){ //custom
		p->dc_level = dc_level;
		p->dc_region_size = (int)pow(2, p->dc_level+1);
		if(p->dc_region_size<min_size)
			G_fatal_error("The most coarse level size must be between %d and pattern size (%d): %d.",min_size,pattern_size,p->dc_region_size);
	}
	else{ //dynamic - if dc_level==0
		p->dc_region_size=2;
		while(p->dc_region_size<=pattern_size){
			p->dc_region_size*=2;
		}
		p->dc_region_size/=2;
		p->dc_level=(int)(log2(p->dc_region_size)-log2(p->dc_base_size)+1);
	}
	p->dc_num_of_quads=p->dc_region_size/p->dc_base_size;

	/*G_message("Number of decomposition levels: %d",p->dc_level); */

	/* end of parameters */

	/* lowering step enables overlaping */

	step=p->dc_region_size;

	/* round to value divided by REGION_SIZE */

	nrows=f->rows - p->dc_region_size + 1;
	ncols=f->cols - p->dc_region_size + 1;

	for(r=1;r<nrows;r+=step)
		for(c=1;c<ncols;c+=step) {
			nulls+=full_decompose_area(l,f,p,r,c,signature);
			count++;
		}

	/* set number of nulls correctly */

	total_cells=count*p->dc_region_size*p->dc_region_size;
	free(p);
	
	if(nulls*2 > total_cells)

		/* mark as null landscape */

		return 0; 
	else
		return 1;
}

int full_decomposition_len(EZGDAL_LAYER **layers, int num_of_layers, ...) {
	int ncats = layers[0]->stats->map_max_val+1;
	
	int dc_level;
	va_list arg_list;
	va_start(arg_list,num_of_layers);
	dc_level = va_arg(arg_list, int);
	va_end(arg_list);
	
	/* # of levels x # of map categories x # of size divisions (hard-coded 3) */
	return dc_level*ncats*3;
}

