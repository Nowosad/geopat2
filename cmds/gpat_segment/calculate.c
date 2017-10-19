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

int compare_by_dist (const void *a, const void *b) {
	return (int)(1000000*(*(DIST **)a)->distance - 1000000*(*(DIST **)b)->distance);
}

int compare_by_dist_desc (const void *a, const void *b) {
	return (int)(1000000*(*(DIST **)b)->distance - 1000000*(*(DIST **)a)->distance);
}

int sort_desc (const void * a, const void * b)
{
   return ( *(int*)b - *(int*)a );
   /* sort descending */
 }

int sort_asc (const void * a, const void * b)
{
   /* sort ascending */
   return ( *(int*)a - *(int*)b );
}

int sort_double_desc (const void * a, const void * b)
{
   /* sort descending */
   return (int)( 1000000*(*(double*)b) - 1000000*(*(double*)a) );
}

int sort_double_asc (const void * a, const void * b)
{
   /* sort ascending */
   return (int)( 1000000*(*(double*)a) - 1000000*(*(double*)b) );
}


double* new_histogram(DATAINFO* d, LOCAL_PARAMS* pars, unsigned index)
{
/* read histogram from data and change data from integer (storage) to double (calculate) */

	static double* ndh; /* new double histogram */
	double div_hist;
	int i,p;
	HISTOGRAM* h=d->histograms;
	long size=(d->size_of_histogram+2)*sizeof(int);
	long offset=size*index;

	fseek(d->fd,offset,SEEK_SET);
	p=(int)fread(d->buffer,size,1,d->fd);

	if(p!=1)
		G_fatal_error(_("Reading histogram data file <%s> error"),d->data_name);

	if((*h->num_of_nulls)==(*h->num_of_features)) {
		ndh=NULL;

		return ndh;
	}

	if((double)(*h->num_of_nulls)/(*h->num_of_features)<pars->null_threshold) {
		div_hist=(double)(*h->num_of_features-*h->num_of_nulls);
		div_hist=1./div_hist;
		ndh=malloc(d->size_of_histogram*sizeof(double));
		for(i=0;i<d->size_of_histogram;++i)
			ndh[i]=div_hist*h->histogram[i];
	}
	else {
		ndh=NULL;
	}

	return ndh;
}

int read_histograms_to_memory(DATAINFO* d, LOCAL_PARAMS* p)
{
	/* for sake of medoids and homogenity it would be better to store it in one place than scattered among areas */
	int i;
	int nrows=d->cell_hd.rows;
	int ncols=d->cell_hd.cols;
	int ncells=nrows*ncols;

	d->all_histograms=malloc(ncells*sizeof(double*));

	/* new histogram returns null if doesn't contain enough data */
	for(i=0;i<ncells;++i)
		d->all_histograms[i]=new_histogram(d, p, i);

	return 0;
	/* for low memory systems it can be easily extended into hdd/ssd swaping system */

}

/* ==================================================================== */

int remove_area(struct area** a)
{
	if((*a)->histogram!=NULL) {
		free((*a)->histogram);
		(*a)->histogram=NULL;
	}
	if((*a)->neighbors!=NULL) {
		clear_queue((*a)->neighbors);
		free((*a)->neighbors);
		(*a)->neighbors=NULL;
	}
	if((*a)->bbox!=NULL) {
		free((*a)->bbox);
		(*a)->bbox=NULL;
	}

	free(*a);
	*a=NULL;

	return 0;
}

int* sample_histogram_ids(struct fifo* queue , int num_of_samples) {

	int i;
	int* ids;
	int* samples=NULL;

	if(num_of_samples<0) /*  to avoid resampling */
		num_of_samples=queue->length;

	if(num_of_samples>queue->length)
		G_fatal_error("Queue of IDS for area %d is too short",queue->first->index);

	ids=convert_queue_to_table(queue);

	if(num_of_samples<queue->length) {
		samples=create_random_sequence(num_of_samples,queue->length);
		for(i=0;i<num_of_samples;++i)
			samples[i]=ids[samples[i]];
		free(ids);
	} else {
		samples=ids;
	}

	return samples;
}

double* use_histogram(DATAINFO* d, int index)
{
	return (index<0)?NULL:d->all_histograms[index];
	/* this function will be extended in the future to get data from SSD or from memory*/
}

int add_histograms(DATAINFO* d, double* o, double* h, int num_of_areas)
{
	double o_div=(num_of_areas-1.)/(double)num_of_areas;
	double h_div=1./(double)num_of_areas;
	int i;

	for(i=0;i<d->size_of_histogram;++i)
		o[i]=o[i]*o_div+h[i]*h_div;

	return 0;
}

double calculate2(HEXGRID* hx, LOCAL_PARAMS* p, double** pair)
{
	/* this function is called whenever previously calculate was called */
	/* we use geometric mean so we want to have distance 1 (maximim) if any of
	 * component is 1
	 * this function calculates similarity (1-dist)
	 * if similarity is 0 (distance is 1) than entire similarity is 0  */
	double result=0, c;
	int i,length=0;
	distance_func* calc=p->calculate;
	int dims[1];

	dims[0]=p->parameters->size_of_histogram;
	if(hx->num_of_subhistograms==1)
		result=1.-calc(pair,2,p->parameters->size_of_histogram,1,dims);
	else {
		int set_0=0;
		for(i=0;i<hx->num_of_subhistograms;++i) {

			double* subpair[]={pair[0]+length,pair[1]+length,NULL};
			p->parameters->size_of_histogram=hx->sh_size_of_histogram[i];

			dims[0]=p->parameters->size_of_histogram;
			c=1.-calc(subpair,2,p->parameters->size_of_histogram,1,dims);
			//c=1.-calc(subpair,p->parameters);
			if(p->all_layers) {
				result=MAX(result,c);
			}
			else {
				set_0=(c<=0)?1:set_0;
				result+=((c<=0)?0:(log(c)*hx->sh_weights[i]));
			}
			length+=hx->sh_size_of_histogram[i];
		}
		if(!p->all_layers)
			result=(set_0)?1:exp(result);
	}

	return 1.-result;
}

int compare_grids_datainfo(DATAINFO* a, DATAINFO* b)
{
  if((a->cell_hd.rows != b->cell_hd.rows) || (a->cell_hd.cols != b->cell_hd.cols))
			G_fatal_error(_("Size of histogram's region in metadata don't match: \n<%s>: %dx%d <%s>: %dx%d"),
										a->data_name,a->cell_hd.cols,a->cell_hd.rows,b->data_name,b->cell_hd.cols,b->cell_hd.rows);
   if((a->cell_hd.ew_res != b->cell_hd.ew_res) || (a->cell_hd.ns_res != b->cell_hd.ns_res))
			G_fatal_error(_("Cell's resolution in metadata don't match: \n<%s>: %fx%f <%s>: %fx%f"),a->data_name,a->cell_hd.ew_res,a->cell_hd.ns_res,b->data_name,b->cell_hd.ew_res,b->cell_hd.ns_res);
   if(a->cell_hd.north != b->cell_hd.north)
			G_fatal_error(_("Map north extent in metadata don't match: \n<%s>: %f <%s>: %f"),a->data_name,a->cell_hd.north,b->data_name,b->cell_hd.north);
   if(a->cell_hd.south != b->cell_hd.south)
			G_fatal_error(_("Map south extent in metadata don't match: \n<%s>: %f <%s>: %f"),a->data_name,a->cell_hd.south,b->data_name,b->cell_hd.south);
   if(a->cell_hd.west != b->cell_hd.west)
			G_fatal_error(_("Map west extent in metadata don't match: \n<%s>: %f <%s>: %f"),a->data_name,a->cell_hd.west,b->data_name,b->cell_hd.west);
   if(a->cell_hd.east != b->cell_hd.east)
			G_fatal_error(_("Map east extent in metadata don't match: \n<%s>: %f <%s>: %f"),a->data_name,a->cell_hd.east,b->data_name,b->cell_hd.east);

  return 0;
}

int get_num_of_grids(char* files[])
{
	FILE* tmp_fd;
	char* tmp_name;
	int i;

	for(i=0;files[i];++i) {
		tmp_name=check_input_names(files[i],"grd");
		tmp_fd=fopen(tmp_name,"r");
			if(tmp_fd==NULL)
		G_fatal_error("Cannot find data file <%s>",tmp_name);
		fclose(tmp_fd);
		free(tmp_name);
		tmp_name=NULL;
		tmp_name=check_input_names(files[i],"hr");
		tmp_fd=fopen(tmp_name,"r");
			if(tmp_fd==NULL)
		G_fatal_error("Cannot find header file <%s>",tmp_name);
		fclose(tmp_fd);
		free(tmp_name);
		tmp_name=NULL;
	}

	return i;
}




