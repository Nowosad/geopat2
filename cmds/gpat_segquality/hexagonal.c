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
#include "local_proto.h"

static int xnextr[4] = { -1, 0, 1,  0 };
static int xnextc[4] = {  0,-1, 0,  1 };

int in_grid(int base, int row_offset, int col_offset, int nrows, int ncols)
{
	int r=base/ncols;
	int c=base%ncols;

	r=r+row_offset;
	if(r<0 || r>=nrows)
		return -1;
	int min_index=r*ncols;
	int max_index=min_index+ncols-1;
	c=c+col_offset;
	int index=r*ncols+c;

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
	int nrows=d->cell_hd.rows;
	int ncols=d->cell_hd.cols;
	int* n=malloc(5*sizeof(int));
	memset(n,-1,5*sizeof(int));

	int r,c,i,j=0;
	r=index/ncols;
	c=index%ncols;
	for(i=0;i<4;++i) {
		if((r+xnextr[i] < 0 || r+xnextr[i] > (nrows-1) || c+xnextc[i] < 0 || c+xnextc[i] > (ncols-1)))
			continue;
		n[j++]=INDEX(r+xnextr[i],c+xnextc[i]);
	}
	qsort(n, 4, sizeof(int), sort_desc);
	return n;
}


int* hex_init_histogram_list(DATAINFO* d,  HEXGRID* hx, int hindex)
{
	int nrows=d->cell_hd.rows;
	int ncols=d->cell_hd.cols;
	int h_r=hindex/hx->ncols;
	int h_c=hindex%hx->ncols;
	int base_index,base_row,max_index,min_index,index,shift=0;
	int i=0,k=0;
	int* list=malloc(5*sizeof(int));

	/* check grid overlap */
	base_index = h_r *2 * ncols + h_c *2; /* base cell, upper left corner */
	base_row=h_r*2;
	shift=h_r%2;

	for(i=0;i<4;++i) {
		min_index=(base_row+(i/2))*ncols;
		max_index=min_index+ncols-1;
		list[i]=-1;

		if(min_index>nrows*ncols-1)
			continue;
		index=base_index+((i/2)*ncols)+(i%2)-shift;

		if((index>=min_index && index<=max_index) && use_histogram(d,index)) {
			list[k++]=index;
//			list[i]=index;
		}
	}
	list[4]=-1; /* close the list */
	if(k==0) {
		free(list);
		list=NULL;
	}
		return list;
}

double* hex_create_histogram(DATAINFO* d, HEXGRID* hx, int index)
{
	int k;
	if(hx->histogram_ids[index]==NULL)
		return NULL;

	double* h;
	int* histogram_ids=hx->histogram_ids[index];
	double* histogram=calloc(d->size_of_histogram,sizeof(double));

	k=0;
	while(histogram_ids[k]>-1) {
		h=use_histogram(d,histogram_ids[k]);
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



struct area* hex_new_area(HEXGRID* hx, unsigned index)
{
	struct area* area=NULL;

	area=malloc(sizeof(struct area));
	area->id=index;
	area->num_of_areas=1;
	area->num_of_samples=1;
	area->num_of_neighbors=0;
	area->done=0;
	area->heterogeneity=0;
	area->isolation=0;

/*	area->hex_ids=malloc(sizeof(int));
	area->hex_ids[0]=index;
	area->sample_ids=area->hex_ids;
	area->clique=area->hex_ids;
*/
	area->hex_ids=NULL;
	area->sample_ids=NULL;

	area->histogram=NULL;
	area->neighbors=NULL;
	area->edges=NULL;
	return area;
}


int hex_remove_area(struct area** a)
{
	if((*a)->histogram!=NULL) {
		free((*a)->histogram);
		(*a)->histogram=NULL;
	}
	if((*a)->neighbors!=NULL) {
		free((*a)->neighbors);
		(*a)->neighbors=NULL;
	}
	if((*a)->edges!=NULL) {
		free((*a)->edges);
		(*a)->edges=NULL;
	}
	if((*a)->hex_ids!=NULL) {
		free((*a)->hex_ids);
		(*a)->hex_ids=NULL;
	}
	if((*a)->sample_ids!=NULL) {
		free((*a)->sample_ids);
		(*a)->sample_ids=NULL;
	}
	return 0;
}

int remove_areas(struct area** areas, int num_of_areas)
{
	int i;
	for(i=0;i<num_of_areas;++i)
		if(areas[i])
			hex_remove_area(areas+i);
	return 0;
}


/* =================================================================== */

HEXGRID* hex_build_topology(DATAINFO* d, LOCAL_PARAMS* p)
{
	int nrows=d->cell_hd.rows;
	int ncols=d->cell_hd.cols;
	int rf=p->quad_mode?1:2;
	int n,i;

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

	hx->nareas=hx->nrows*hx->ncols;
	hx->reduction=rf;
	hx->size_of_supermotifel=rf*rf;
	hx->size_of_histogram=d->size_of_histogram;
	hx->histogram_ids=malloc(hx->nareas*sizeof(int*));
	hx->histograms=malloc(hx->nareas*sizeof(double*));
	hx->hex_neigborhoods=malloc(hx->nareas*sizeof(int*));

	if(p->quad_mode) {
		for(i=0;i<hx->nareas;++i) {
			hx->histograms[i]=NULL;
			hx->hex_neigborhoods[i]=NULL;
			hx->histogram_ids[i]=NULL;

			hx->histograms[i]=use_histogram(d,i);
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
			int* list=hex_init_histogram_list(d,hx,i);
			if(list) {
				hx->histogram_ids[i]=list;
				hx->histograms[i]=hex_create_histogram(d,hx,i);
			}
		}
	}

	if(p->quad_mode){
		for(i=0;i<hx->nareas;++i)
			hx->hex_neigborhoods[i]=hex_get_quad_neighborhood(d,hx,i);
	} else {
		for(i=0;i<hx->nareas;++i)
			hx->hex_neigborhoods[i]=hex_get_neighborhood(d,hx,i);
	}
	return hx;
}


int hex_remove_hexgrid(HEXGRID* hx)
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

	return 0;
}

