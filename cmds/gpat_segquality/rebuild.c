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

int sample_area(HEXGRID* hx, LOCAL_PARAMS* pars, struct area* a)
{
	int i;
	int* samples=NULL;

	if(a->num_of_areas<=(int)pars->sampling_threshold*1.1||pars->sampling_threshold==0)
		a->num_of_samples=a->num_of_areas;
	else
		a->num_of_samples=pars->sampling_threshold;

	a->sample_ids=malloc(a->num_of_samples*sizeof(int));
	if(a->num_of_samples<a->num_of_areas) {
		samples=create_random_sequence(a->num_of_samples,a->num_of_areas);
		for(i=0;i<a->num_of_samples;++i) {
			a->sample_ids[i]=a->hex_ids[samples[i]];
		}
		free(samples);

	} else {
		for(i=0;i<a->num_of_samples;++i) /* just copy. Think about pointer */
			a->sample_ids[i]=a->hex_ids[i];
	}
	return 0;
}

struct area** build_areas(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* p, int* results)
{
	int i,index=0;
	struct area** areas=NULL;
	areas=(struct area**)calloc(hx->nareas,sizeof(struct area*));
	struct fifo** hex_ids=calloc(hx->nareas,sizeof(struct fifo*));

	/* fill new areas */
	for(i=0;i<hx->nareas;++i) {
		if(!Rast_is_c_null_value(results+i) && results[i]>-1) {
			index=results[i];
			if(areas[index]==NULL) {
				areas[index]=hex_new_area(hx,index);
				areas[index]->id=index;
				hex_ids[index]=init_queue(i);
			} else {
				add_node(hex_ids[index],i);
			}
		}
	}
	for(i=0;i<hx->nareas;++i) {
		if(areas[i]) {
			areas[i]->hex_ids=convert_queue_to_table(hex_ids[i]);
			areas[i]->num_of_areas=queue_length(hex_ids[i]);
			clear_queue(hex_ids[i]);
			sample_area(hx,p,areas[i]);
		}
	}
	free(hex_ids);
	return areas;
}

int build_neighbors_queue(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results)
{
	/* counts number of possible connections and put it to the queues
	 * to reduce number of calculation before sorting we only join segments with smaller index
	 * to segment with higher index
	 * It allow to avoid double calculations of simetrical pairs: for example 12--13 is equal to pair 13--12
	 * we use only the first so first condition is: index_segment<target_segment */

	 /* important note: length of the border is stored in queue */

	int i,index,target,index_segment,target_segment;
	struct fifo** neighbors=calloc(hx->nareas,sizeof(struct fifo*));

	for(index=0;index<hx->nareas;++index) {

		if(Rast_is_c_null_value(results+index))
			continue;

		index_segment=results[index];
		if(areas[index_segment]==NULL)
			continue;

		if(!neighbors[index_segment])
			neighbors[index_segment]=init_queue(index_segment);

		target_segment=-1;
		i=0;

		while(hx->hex_neigborhoods[index][i] >-1) {
			target=hx->hex_neigborhoods[index][i];
			i++;
			target_segment=results[target];
			if(Rast_is_c_null_value(&target_segment))
				continue;
			if(areas[target_segment]==NULL) /* added ealier but yet not used */
				G_fatal_error("AREA does not exist"); //MOZE ZBEDNE


			if(index_segment!=target_segment)
				if(!in_queue(neighbors[index_segment],target_segment)) {
					add_node(neighbors[index_segment],target_segment);
				}
		}
	}
	for(index=0;index<hx->nareas;++index) {
		if(areas[index]) {
			areas[index]->neighbors=convert_queue_to_table(neighbors[index]);
			areas[index]->edges=convert_queue_to_edges(neighbors[index]);
			areas[index]->num_of_neighbors=queue_length(neighbors[index]); /* is one up */
			clear_queue(neighbors[index]);
		}
	}
	free(neighbors);
	return 0;
}

int* create_segment_map(HEXGRID* hx, LOCAL_PARAMS* p, EZGDAL_LAYER* l)
{
/* creates right map of segments depending on whether it's brick or rook topology */
	int r,c,i,j;
	int nrows = l->rows;
	int ncols = l->cols;
	int segment,fragment;
	int* input_data = malloc(nrows*ncols*sizeof(int));
	int* segments;

	/* read input raster to memory */
	for(r=0; r<nrows; r++) {
		ezgdal_read_buffer(l, r);
		for(c=0; c<ncols; c++)
			if(ezgdal_is_null(l, l->buffer[c]))
				input_data[r*ncols+c]=0x80000000;
			else
				input_data[r*ncols+c] = (int)(l->buffer[c]);
	}

	if(p->quad_mode)
		return input_data;

	segments = malloc(hx->nareas*sizeof(int));
	Rast_set_c_null_value(segments,hx->nareas);

	for(i=0;i<hx->nareas;++i) {
		j=0;
		if(hx->histogram_ids[i]) {
			segment=input_data[hx->histogram_ids[i][j]];
			while(hx->histogram_ids[i][j]>-1) {
				fragment=input_data[hx->histogram_ids[i][j]];
				if(!Rast_is_c_null_value(&segment) && !Rast_is_c_null_value(&fragment) && segment!=fragment)
					G_fatal_error(_("Segmentation map has no brick topology: %d %d"),segment,input_data[hx->histogram_ids[i][j]]);
				j++;
			}
			segments[i]=segment;
		}
	}

	free(input_data);
	return segments;
}

