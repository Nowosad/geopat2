/****************************************************************************
 *
 * MODULE:	segmentation of motifels grid
 * AUTHOR(S):	Jaroslaw Jasiewicz, Jacek Niesterowicz, Tomasz Stepinski
 * PURPOSE:	information retrival using categorical maps:
 *		compares grid of histograms
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include "local_proto.h"

double recalculate_threshold(HEXGRID* hx, LOCAL_PARAMS* pars, struct area* a)
{
	int i,j,num_of_a_samples;
	int* a_samples;
	int num_of_pairs=0;
	double distance=0;
	double mean_distance=0;
	double stddev_distance=0;

	/* sampling */
	if(a->num_of_areas<=(int)pars->sampling_threshold*1.1||pars->sampling_threshold==0)
		num_of_a_samples=a->num_of_areas;
	else
		num_of_a_samples=pars->sampling_threshold;
	a_samples=sample_histogram_ids(a->hex_ids,num_of_a_samples);

	if(num_of_a_samples==1) {
		a->homogenity=0;

		return 0;
	}

	for(i=0;i<num_of_a_samples;++i) {
		for(j=(i+1);j<num_of_a_samples;++j) {
			double* pair[]={NULL,NULL,NULL};
			pair[0]=hex_use_histogram(hx,a_samples[i]);
			pair[1]=hex_use_histogram(hx,a_samples[j]);
			distance=calculate2(hx,pars,pair);
			mean_distance+=distance;
			stddev_distance+=(distance*distance);
			num_of_pairs++;
		}
	}
	if(num_of_pairs>0) {
		mean_distance/=(double)num_of_pairs;
		stddev_distance/=(double)num_of_pairs;
		stddev_distance=sqrt(stddev_distance-(mean_distance*mean_distance));
	}

	free(a_samples);

	return mean_distance+stddev_distance;
	//return mean_distance;
}


int add_nodes_to_queue(int* neighbors, int* results, struct fifo* queue)
{
	/* adds all areas surrounding given area to queue */
	int i=0,target=0;

	while(neighbors[i]>-1) {
		target=neighbors[i];
		i++;
		if(results[target]>-1||Rast_is_c_null_value(results+target))
			continue;
		if(!in_queue(queue,target))
			add_node(queue,target);
	}

	return 0;
}

int find_most_similar(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, struct fifo* queue, int index)
{
	/* we are inside queue, ut uses queue api */
	struct node* next;
	double min_distance=areas[index]->similarity_threshold;
	double distance;
	int min_index=-1;

	next=queue->first;

	if(areas[index]->num_of_areas<1.1*pars->sampling_threshold||pars->sampling_threshold==0) {
		while(next) {
			distance=hex_calculate_linkage(hx,pars,areas[index],areas[next->index]);
			if(distance<min_distance) {
				min_distance=distance;
				min_index=next->index;
			}
			next=next->next_node;
		}
		if(min_index>-1)
			find_and_remove_node(queue,min_index);
		return min_index;
	} else {
		while(next) {
			distance=hex_calculate_linkage(hx,pars,areas[index],areas[next->index]);
			if(distance<areas[index]->similarity_threshold) {
				min_index=next->index;
				remove_node(queue);
				return min_index;
			}
			next=next->next_node;
			remove_node(queue);
		}
	}

	return min_index;
}

int build_segment(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results, int segment)
{
	int target;
	struct fifo* queue=init_queue(segment);

	add_nodes_to_queue(hx->hex_neigborhoods[segment],results,queue);
	remove_node(queue); /* remove first element -- a segment */
	areas[segment]->modified=0;
	while (queue->first) {
		if(!areas[segment]->modified && areas[segment]->num_of_areas>1.1*pars->sampling_threshold && pars->sampling_threshold>0) {
			double segment_heterogeneity=recalculate_threshold(hx,pars,areas[segment]);
			if(segment_heterogeneity<areas[segment]->similarity_threshold)
				areas[segment]->similarity_threshold=segment_heterogeneity;
				//areas[segment]->similarity_threshold=(segment_heterogeneity+areas[segment]->similarity_threshold)/2;
			areas[segment]->modified=1;
		}
		target=find_most_similar(hx,pars,areas,queue,segment);
			if(target>-1) {
				results[target]=segment;
				if(queue->first)
					add_nodes_to_queue(hx->hex_neigborhoods[target],results,queue);
				hex_join_areas(areas[segment],areas[target]);
				areas[target]=NULL;
			} else
				clear_queue(queue); /* none of element is below the threshold */
	}
	if(queue)
		free(queue); /* queue is to be empty maybe check before free? */
	areas[segment]->modified=0;

	return 0;
}

int hex_region_growing(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results, unsigned* seeds, unsigned num_of_seeds)
{
	unsigned int i,num_of_segments=0;
	int segment;

	G_message("Growing segments...");
	Rast_set_c_null_value(results,hx->nareas);
	for(i=0;i<num_of_seeds;++i)
		results[seeds[i]]=-1;

	for(i=0;i<num_of_seeds;++i){

		segment=seeds[i];
		if(areas[segment]==NULL)
			continue;

		results[segment]=segment;
		num_of_segments++;
		build_segment(hx,pars,areas,results,segment);

	if(num_of_segments<250||num_of_segments%50==0)
		fprintf(stderr, "%08d\b\b\b\b\b\b\b\b", num_of_segments);
	} /* end seeds */
	fprintf(stderr, "\b\b\b\b\b\b\b\b");
	G_message("Number of grown segments: %d", num_of_segments);
	free(seeds);

	return 0;
}


