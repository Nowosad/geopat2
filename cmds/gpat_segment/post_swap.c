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

int swap_cell(struct area** areas, int *results, int index, int segment, int target)
{
	find_and_remove_node(areas[segment]->hex_ids,index);
	add_node(areas[target]->hex_ids,index);
	areas[target]->num_of_areas++;
	areas[segment]->num_of_areas--;
	results[index]=target;
	if(areas[segment]->num_of_areas==0) {
		remove_area(areas+segment);
		areas[segment]=NULL;
		return 1;
	}
	return 0;
}

int* is_border(HEXGRID* hx, LOCAL_PARAMS* p, struct area** areas, int *results, int index)
{
	int i,j,k;
	int found,target_segment,target;
	int* neighbors=malloc(7*sizeof(int));
	memset(neighbors,-1,7*sizeof(int));
	k=0;

	i=0;
	while(hx->hex_neigborhoods[index][i] >-1) {
		target=hx->hex_neigborhoods[index][i];
		target_segment=results[target];
		if(!Rast_is_c_null_value(results+target) && results[index]!=target_segment) {
			found=0;
			for(j=0;j<=k;++j) /* check if not already added */
				if(neighbors[j]==target_segment) {
					found=1;
					break;
				}
			if(!found)
				neighbors[k++]=target_segment;
		}
	i++;
	}
	if(k==0) {
		free(neighbors);
		return NULL;
	}
	return neighbors;
}

int do_swap_cells(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results)
{

	int i,target_segment,index_segment,index,cur_segment;
	double segment_distance,target_distance,distance;
	struct area* current_area=NULL;
	int* neighbors=NULL;
	int num_of_swapped_cells=0;

	double swapping_threshold=pars->lower_similarity_threshold/10;
	for(index=0;index<hx->nareas;++index) {
		index_segment=results[index];
		if(Rast_is_c_null_value(results+index))
			continue;
		target_distance=2;
		target_segment=-1;
		neighbors=is_border(hx,pars,areas,results,index);
		if(!neighbors)
			continue;
		current_area=hex_new_area(hx,index);

		if(queue_length(areas[index_segment]->hex_ids)==1)
			continue;
		else {
			find_and_remove_node(areas[index_segment]->hex_ids,index); /* we need to remove tested motifel from segment to avoid self-similarity */
			areas[index_segment]->num_of_areas--;
			segment_distance=hex_calculate_linkage(hx,pars,areas[index_segment],current_area);
			add_node(areas[index_segment]->hex_ids,index); /* and we add it back */
			areas[index_segment]->num_of_areas++;
		}
		if(segment_distance<swapping_threshold) { /* are originally too similar to be swapped */
			free(neighbors);
			remove_area(&current_area);
			continue;
		}
		i=0;
		target_distance=2;
		target_segment=-1;
		while(neighbors[i]>-1) {
			cur_segment=neighbors[i];
			distance=hex_calculate_linkage(hx,pars,areas[cur_segment],current_area);
			if(distance<target_distance) {
				target_distance=distance;
				target_segment=cur_segment;
			}
			i++;
		}

		if(target_segment>=0) {
			if((segment_distance-target_distance)>swapping_threshold) {
				swap_cell(areas,results,index,index_segment,target_segment);
				num_of_swapped_cells++;
			}
		}
		if(current_area) {
			remove_area(&current_area);
			current_area=NULL;
		}
	}
	fprintf(stderr, "%08d\b\b\b\b\b\b\b\b", num_of_swapped_cells);
	return num_of_swapped_cells;
}

int add_nodes_to_queue_id(HEXGRID* hx, int* results, struct fifo* queue, unsigned index)
{
	/* adds all areas surrounding given areaa to queue */

	int i=0, target;
	while(hx->hex_neigborhoods[index][i] >-1) {
		target=hx->hex_neigborhoods[index][i];
		if(results[target]==results[index] && !in_queue(queue,target))
			add_node(queue,target);
		i++;
	}
	return 0;
}

struct area* rebuild_area(HEXGRID* hx, LOCAL_PARAMS* pars, int* results, unsigned index)
{
	struct area* area=hex_new_area(hx,index); /* new area inits queue with index as first element */;
	struct node* current_node=area->hex_ids->first;

	while(current_node) {
		add_nodes_to_queue_id(hx,results,area->hex_ids,current_node->index);
		current_node=current_node->next_node;
	}
	area->num_of_areas=queue_length(area->hex_ids);
	return area;
}

int rebuild_areas(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int* results)
{
	/* clear areas */
	int i,j,k=1;
	int* list_of_cells=NULL;

	fprintf(stderr, "...Re-build segments...\n");
	for(i=0;i<hx->nareas;++i) {
		if(areas[i]) {
			remove_area(&(areas[i]));
		}
		areas[i]=NULL;
	}

	/* build new areas */
	for(i=0;i<hx->nareas;++i)
		if(!Rast_is_c_null_value(results+i) && results[i]>-1) {
			if(areas[i]==NULL) {
				areas[i]=rebuild_area(hx,pars,results,i);
				list_of_cells=convert_queue_to_table(areas[i]->hex_ids);
				for(j=0;j<areas[i]->num_of_areas;++j) {
					results[list_of_cells[j]]=-(i+1); /* to avoid -0 */ //re-numeracja do poprawy nie jest potrzebna
				}
				free(list_of_cells);
				list_of_cells=NULL;
				areas[i]->id=k++;
			}
		}

	/* re-numerate results layer */
	for(i=0;i<hx->nareas;++i)
		if(!Rast_is_c_null_value(results+i) && results[i]<0)
			results[i]=-(results[i]+1);

	return 0;
}

int swap_areas(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results)
{
	int total_num_of_swapped_cells=0;
	int num_of_swapped_cells;

	fprintf(stderr,"...Swap segments...\n");
	do{
		num_of_swapped_cells=0;
		num_of_swapped_cells=do_swap_cells(hx,pars,areas,results);
		total_num_of_swapped_cells+=num_of_swapped_cells;
	} while(num_of_swapped_cells>(pars->swap_threshold*hx->nareas)); /* free parameter  move to menu instead of geometry*/

	fprintf(stderr, "\b\b\b\b\b\b\b\b");
	fprintf(stderr,"Num of swapped areas: %d... Done\n", total_num_of_swapped_cells);
	rebuild_areas(hx,pars,areas,results);
	return 0;
}
