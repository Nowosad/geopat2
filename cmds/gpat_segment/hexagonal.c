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
#include "local_proto.h"

static int xnextr[4] = { -1, 0, 1,  0 };
static int xnextc[4] = {  0,-1, 0,  1 };

int in_grid(int base, int row_offset, int col_offset, int nrows, int ncols)
{
	int min_index, max_index, index;
	int r=base/ncols;
	int c=base%ncols;

	r=r+row_offset;
	if(r<0 || r>=nrows)
		return -1;
	min_index=r*ncols;
	max_index=min_index+ncols-1;
	c=c+col_offset;
	index=r*ncols+c;

	if(index<min_index || index>max_index)
		return -1;
	return index;
}

double* hex_use_histogram(HEXGRID* hx, int index)
{
	return hx->histograms[index];
}

int* hex_get_neighborhood(DATAINFO* d, HEXGRID* hx, int index)
{
	/* this function returns list of indexes for direct neighborhood in pseudo-hex topology grid */

	int* n = malloc(7*sizeof(int));
	int h_row = index / hx->ncols;
	int i;
	int shift = h_row % 2; /* shift 0 for odd and -1 for even row */

/*      0 | 1
	 * 2 | x | 3
	 *   4 | 5     */

	n[0]=in_grid(index,-1,-shift,hx->nrows,hx->ncols);
	n[1]=in_grid(index,-1,1-shift,hx->nrows,hx->ncols);
	n[2]=in_grid(index,0,-1,hx->nrows,hx->ncols);
	n[3]=in_grid(index,0,1,hx->nrows,hx->ncols);
	n[4]=in_grid(index,1,-shift,hx->nrows,hx->ncols);
	n[5]=in_grid(index,1,1-shift,hx->nrows,hx->ncols);
	n[6]=-1;

	for(i=0;i<6;++i) {
		if(n[i]>-1 && hx->histograms[n[i]]==NULL)
			n[i]=-1;
	}
	qsort(n, 6, sizeof(int), sort_desc);

	return n; /* to AREA */
}

int* hex_get_quad_neighborhood(DATAINFO* d, HEXGRID* hx, int index)
{
	/* creates neighborhood for quad (rook) topology */
	int r,c,i,j=0;
	int nrows=d->cell_hd.rows;
	int ncols=d->cell_hd.cols;
	int* n=malloc(5*sizeof(int));

	memset(n,-1,5*sizeof(int));

	r=index/ncols;
	c=index%ncols;
	for(i=0;i<4;++i) {
		if((r+xnextr[i] < 0 || r+xnextr[i] > (nrows-1) || c+xnextc[i] < 0 || c+xnextc[i] > (ncols-1)))
			continue;
		if(!use_histogram(d,INDEX(r+xnextr[i],c+xnextc[i])))
			continue;
		n[j++]=INDEX(r+xnextr[i],c+xnextc[i]);
	}
	qsort(n, 4, sizeof(int), sort_desc);

	return n;
}


int* hex_init_histogram_list(DATAINFO** d, LOCAL_PARAMS* p, HEXGRID* hx, int hindex)
{
	int nrows=d[0]->cell_hd.rows;
	int ncols=d[0]->cell_hd.cols;
	int h_r=hindex/hx->ncols;
	int h_c=hindex%hx->ncols;
	int base_index,base_row,max_index,min_index,index,shift=0;
	int i=0,k=0,j,rf=hx->reduction; /* rf - reduction factor */
	int* list=malloc(hx->size_of_supermotifel*sizeof(int));

	/* check grid overlap */
	base_row=h_r*rf;
	shift=(h_r%2)?rf/2:0; /* we shift every second region */
	base_index = base_row * ncols  + h_c * rf; /* base cell, upper left corner */

	for(i=0;i<hx->size_of_supermotifel;++i) {
		min_index=(base_row+(i/rf))*ncols;
		max_index=min_index+ncols;
		list[i]=-1;

		if(min_index>nrows*ncols-1 || max_index<0) { /* out of grid */
			continue;
		}
		index=base_index+((i/rf)*ncols)+(i%rf)-shift;

		if(index<min_index || index>=max_index)
			continue;

		int complete=1;
		for(j=0;j<hx->num_of_subhistograms;++j) {
			if(!use_histogram(d[j],index))
				complete=0;
		}
		if(complete) {
			list[i]=index;
			k++;
		}

	}
	if(k==0) {
		free(list);
		list=NULL;
	}
	return list;
}

double* hex_create_histogram(DATAINFO* d, HEXGRID* hx, int index)
{
	int i,k=0;
	double* h;
	int* histogram_ids;
	double* histogram;

	if(hx->histogram_ids[index]==NULL)
		return NULL;

	histogram_ids=hx->histogram_ids[index];
	histogram=calloc(d->size_of_histogram,sizeof(double));

	for(i=0;i<hx->size_of_supermotifel;++i) {
		h=use_histogram(d,histogram_ids[i]);
		if(h) {
			k++;
			add_histograms(d,histogram,h,k);
		}
	}
	if(!k) {
		free(histogram);
		return NULL;
	}

	return histogram;
}

double* hex_join_histograms(DATAINFO** d, HEXGRID* hx, int index)
{
	/* to unified calculations we create one histogram which is passed to calculate function
	 * elements of histograms can be accessed by slicing
	 */

	int i;
	double* histogram;
	double* h;
	int start=0;

	if(hx->histogram_ids[index]==NULL)
		return NULL;

	if(hx->num_of_subhistograms==1) {
		histogram=hex_create_histogram(d[0],hx,index);
		return histogram;
	}

	histogram=malloc(hx->size_of_histogram*sizeof(double));
	for(i=0;i<hx->num_of_subhistograms;++i) {
		h=hex_create_histogram(d[i],hx,index);
		memcpy(histogram+start,h,hx->sh_size_of_histogram[i]*sizeof(double));
		start+=hx->sh_size_of_histogram[i];
		free(h);
	}

	return histogram;
}

double* hex_join_quad_histograms(DATAINFO** d, HEXGRID* hx, int index)
{
	/* difference between quad and hex: hex requires creation of new histogram
	 * while quad do not so we use function use_histogram
	 */
	double* h;
	int i, start=0;
	double* histogram;

	for(i=0;i<hx->num_of_subhistograms;++i) {
		h=use_histogram(d[i],index);
		if(h==NULL)
			return NULL;
	}

	histogram=malloc(hx->size_of_histogram*sizeof(double));
	for(i=0;i<hx->num_of_subhistograms;++i) {
		h=use_histogram(d[i],index);
		memcpy(histogram+start,h,hx->sh_size_of_histogram[i]*sizeof(double));
		start+=hx->sh_size_of_histogram[i];
	}

	return histogram;
}

struct area* hex_new_area(HEXGRID* hx, unsigned index)
{
	struct area* area=NULL;

	area=malloc(sizeof(struct area));
	area->id=index;
	area->hex_ids=init_queue(index);
	area->num_of_areas=1;
	area->edge=0;
	area->modified=1;
	area->homogenity=0;

	area->histogram=NULL;
	area->bbox=NULL;
	area->neighbors=NULL;

	return area;
}


/* =================================================================== */

HEXGRID* hex_build_topology(DATAINFO** d, LOCAL_PARAMS* p, int num_of_layers, char* weights)
{
	int nrows=d[0]->cell_hd.rows;
	int ncols=d[0]->cell_hd.cols;
	int rf=p->quad_mode?1:(int)pow(2.,p->reduction+1);
	int n,i,size_of_histogram=0;

	HEXGRID* hx=malloc(sizeof(HEXGRID));

	if(p->quad_mode) {
		hx->nrows = nrows;
		hx->ncols = ncols;
	} else {
		n=nrows/rf;
		hx->nrows = n*rf<nrows?n+2:n+1;
		n=ncols/rf;
		hx->ncols = n*rf<ncols?n+2:n+1;
		hx->level=-1;
	}

	hx->sh_size_of_histogram=malloc(num_of_layers*sizeof(int));
	for(i=0;i<num_of_layers;++i) {
		hx->sh_size_of_histogram[i]=d[i]->size_of_histogram;
		size_of_histogram+=d[i]->size_of_histogram;
	}

	hx->sh_weights=parse_weights(num_of_layers,weights);
	hx->num_of_subhistograms=num_of_layers;

	hx->size_of_histogram=size_of_histogram; /* if one histogram only we will use only this value */
	hx->nareas=hx->nrows*hx->ncols;
	hx->reduction=rf;
	hx->size_of_supermotifel=rf*rf;
	hx->thresholds=malloc(hx->nareas*sizeof(double));
	hx->histogram_ids=malloc(hx->nareas*sizeof(int*));
	hx->histograms=malloc(hx->nareas*sizeof(double*));
	hx->hex_neigborhoods=malloc(hx->nareas*sizeof(int*));
	hx->low_res_grid=NULL;

	/* QUAD uses same technology as brick: list of elements but list has only one element */

	if(p->quad_mode) {
		for(i=0;i<hx->nareas;++i) {
			hx->histograms[i]=NULL;
			hx->hex_neigborhoods[i]=NULL;
			hx->histogram_ids[i]=NULL;

			hx->histograms[i]=hex_join_quad_histograms(d,hx, i);
			if(hx->histograms[i]) {
				hx->histogram_ids[i]=malloc(2*sizeof(int));
				hx->histogram_ids[i][0]=i;
				hx->histogram_ids[i][1]=-1;
			}
		}
	} else {
		for(i=0;i<hx->nareas;++i) {
			hx->histograms[i]=NULL;
			hx->hex_neigborhoods[i]=NULL;
			hx->histogram_ids[i]=NULL;
			int* list=hex_init_histogram_list(d,p,hx,i);
			if(list) {
				hx->histogram_ids[i]=list;
				hx->histograms[i]=hex_join_histograms(d,hx,i);
			}
		}
	}

	if(p->quad_mode){
		for(i=0;i<hx->nareas;++i)
			hx->hex_neigborhoods[i]=hex_get_quad_neighborhood(d[0],hx,i);
	} else {
		for(i=0;i<hx->nareas;++i)
			hx->hex_neigborhoods[i]=hex_get_neighborhood(d[0],hx,i);
	}

	return hx;
}


int* hex_remove_hexgrid(HEXGRID* hx)
{
	int i;

	if(hx->histogram_ids) {
		for(i=0;i<hx->nareas;++i)
			if(hx->histogram_ids[i]) {
				free(hx->histogram_ids[i]);
				hx->histogram_ids[i]=NULL;
			}
		free(hx->histogram_ids);
		hx->histogram_ids=NULL;
	}

	if(hx->histograms) {
		for(i=0;i<hx->nareas;++i)
			if(hx->histograms[i]) {
				free(hx->histograms[i]);
				hx->histograms[i]=NULL;
			}
		free(hx->histograms);
		hx->histograms=NULL;
	}

	if(hx->hex_neigborhoods) {
		for(i=0;i<hx->nareas;++i)
			if(hx->hex_neigborhoods[i]) {
				free(hx->hex_neigborhoods[i]);
				hx->hex_neigborhoods[i]=NULL;
			}
		free(hx->histograms);
		hx->histograms=NULL;
	}

	if(hx->thresholds) {
		free(hx->thresholds);
		hx->thresholds=NULL;
	}

	return hx->low_res_grid;
}

struct area** hex_build_areas(DATAINFO** d, HEXGRID* hx, LOCAL_PARAMS* p)
{
	/* build also histograms and neigborhoods */
	int i;
	struct area** areas;

	areas=malloc(hx->nareas*sizeof(struct area*));

	for(i=0;i<hx->nareas;++i) {
		if(hx->histograms[i])
			areas[i]=hex_new_area(hx,i);
		else
			areas[i]=NULL;
	}

	return areas;
}

int* hex_init_results(HEXGRID* hx)
{

	int i;
	int* results=(int*)malloc(hx->nareas*sizeof(int));

	Rast_set_c_null_value(results,hx->nareas);
	for(i=0;i<hx->nareas;++i) {
		if(hx->histograms[i])
			results[i]=i;
	}

	return results;
}

int hex_reclass(HEXGRID* hx, struct area** areas)
{
	int i;
	int k=1;

	for(i=0;i<hx->nareas;++i) {
		if(areas[i]!=NULL) {
			areas[i]->id=k++; /* reclass */
		}
	}
	G_message("Final number of segments: %d",k-1);

	return 0;
}

int* hex_create_segment_map(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* p, struct area** areas)
{
	int i,j,k;
	int nrows=d->cell_hd.rows;
	int ncols=d->cell_hd.cols;
	int ncells=nrows*ncols;
	int* results=malloc(ncells*sizeof(int));
	struct fifo* areas_ids=NULL;
	int* ids;
	int num_of_histograms;
	int* histogram_ids;

	Rast_set_c_null_value(results,ncells);
	/* we switch back to the quad topology */

	for(i=0;i<hx->nareas;++i) {
		if(areas[i]) {
			ids=convert_queue_to_table(areas[i]->hex_ids);
			num_of_histograms=queue_length(areas[i]->hex_ids);
			for (j=0;j<num_of_histograms;++j) {
				histogram_ids=hx->histogram_ids[ids[j]];
				if(histogram_ids==NULL) /* shoudn't happen */
					G_fatal_error("HEXconvert:NULL IDS");

				for(k=0;k<hx->size_of_supermotifel;++k) {
					if(histogram_ids[k]>-1) {
						results[histogram_ids[k]]=areas[i]->id;
						if(!areas_ids)
							areas_ids=init_queue(histogram_ids[k]);
						else
							add_node(areas_ids,histogram_ids[k]);
					}
				}
			}
			clear_queue(areas[i]->hex_ids);
			areas[i]->hex_ids=areas_ids;
			areas_ids=NULL;
			areas[i]->num_of_areas=queue_length(areas[i]->hex_ids);
			free(ids);
		}
	}
	return results;
}

int hex_build_neighbors_queue(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results, int remove)
{
	/* counts number of possible connections and put it to the queues
	 * to reduce number of calculation before sorting we only join segments with smaller index
	 * to segment with higher index
	 * It allow to avoid double calculations of simetrical pairs: for example 12--13 is equal to pair 13--12
	 * we use only the first so first condition is: index_segment<target_segment if we want all connections than index_segment!=target_segment */

	int i=0,index,target,index_segment,target_segment;

	for(index=0;index<hx->nareas;++index) {

		if(Rast_is_c_null_value(results+index))
			continue;

		index_segment=results[index];
		if(areas[index_segment]==NULL)
			continue;

		if(areas[index_segment]->neighbors==NULL)
			continue;

		target_segment=-1;
		i=0;
		while(hx->hex_neigborhoods[index][i] >-1) {
			target=hx->hex_neigborhoods[index][i];
			i++;
			target_segment=results[target];
			if(areas[target_segment]==NULL) /* added ealier but yet not used */
				G_fatal_error("AREA does not exist"); //MOZE ZBEDNE

			if(remove) {
				if(index_segment<target_segment)
					if(!in_queue(areas[index_segment]->neighbors,target_segment)) {
						add_node(areas[index_segment]->neighbors,target_segment);
						areas[index_segment]->num_of_neighbors++;
					}
			} else {
				if(index_segment!=target_segment)
					if(!in_queue(areas[index_segment]->neighbors,target_segment)) {
						add_node(areas[index_segment]->neighbors,target_segment);
						areas[index_segment]->num_of_neighbors++;
					}
			}
		}
	}

	return 0;
}

int hex_recalculate_neighbors(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int* results)
{
	int i;

	for(i=0;i<hx->nareas;++i) {
		if(areas[i]) {
			if(areas[i]->neighbors!=NULL) {
				clear_queue(areas[i]->neighbors); /* clear all areas which contain queue */
				areas[i]->neighbors=NULL;
				areas[i]->num_of_neighbors=0;
			}
			areas[i]->neighbors=init_queue(i); /* first element of queue is area, next elements will be added after it! */
			areas[i]->done=0;
		}
	}
	hex_build_neighbors_queue(hx,pars,areas,results,0); /* do not remove  dual pairs */

	return 0;
}
