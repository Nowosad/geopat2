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

double m_find_distance(HEXGRID* hx, LOCAL_PARAMS* pars, struct area* index, struct area* target)
{
	double tmp_dist=hex_calculate_linkage(hx,pars,index,target);

	if(tmp_dist>pars->upper_similarity_threshold)
		return 2.; /* eliminate */

	return tmp_dist;
}

int hex_minarea (HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results)
{

	int i,j,num_of_pairs=0;
	int index_segment,target_segment;
	int num_of_removed_areas=0;
	int num_of_minareas=0;
	DIST* dists;
	DIST** pointers_to_dists;
	int previous_num_of_removed_areas=0;
	int x=0;

	G_message("Removing small areas by joining...");
	for(i=0;i<hx->nareas;++i) {
		if(areas[i]!=NULL) {
			if(areas[i]->neighbors!=NULL) {
				clear_queue(areas[i]->neighbors); /* clear all areas which contain queue */
				areas[i]->neighbors=NULL;
				areas[i]->num_of_neighbors=0;
			}
			if(areas[i]->num_of_areas<pars->minarea) {
				areas[i]->neighbors=init_queue(i); /* first element of queue is area, next elements will be added after it! */
				areas[i]->done=0;
				num_of_minareas++;
			}
		}
	}
	if(!num_of_minareas) {
		G_message("No small area to remove");
		return 0;
	}
	fprintf(stderr, "...Buliding graph...");
	/* counts number of possible connections and put it to the queues */
	hex_build_neighbors_queue(hx,pars,areas,results,0);

	for(i=0;i<hx->nareas;++i)
		if(areas[i]!=NULL && areas[i]->neighbors && areas[i]->neighbors->length>1) /* first element is an index, so one element means segments has no joins */
			num_of_pairs+=(areas[i]->neighbors->length-1);

	dists=(DIST*)G_malloc(sizeof(DIST)*num_of_pairs);
	pointers_to_dists=(DIST**)G_malloc(sizeof(DIST*)*num_of_pairs);

	/* convert list into table, maybe not very clear */
	int k=0; /* index of pairs */
	for(i=0;i<hx->nareas;++i) {
		if(areas[i]!=NULL && areas[i]->neighbors && areas[i]->neighbors->length>0) { /* as long as there are some neighbors */
			index_segment=pop_node(areas[i]->neighbors);
			while(areas[index_segment]->neighbors->first) {
				target_segment=pop_node(areas[index_segment]->neighbors); /* we remove all targets */
				dists[k].index=index_segment;
				dists[k].result=target_segment; /* use results as targets */
				dists[k].distance=m_find_distance(hx,pars,areas[index_segment],areas[target_segment]);
				*(pointers_to_dists+k)=dists+k;
				k++;
				if(k>num_of_pairs)
					G_fatal_error("Cannot remove areas, wrong number of pairs");
			}
			free(areas[index_segment]->neighbors);
			areas[index_segment]->neighbors=NULL;
		}
	}

	qsort(pointers_to_dists, num_of_pairs, sizeof(DIST*), compare_by_dist);
	fprintf(stderr, " Done \n");

	while (pointers_to_dists[0]->distance<=1) { /* areas are sorted ASC by distance (!NOT SIMILARITY!!) all others have distance grater than threshold so we can break for loop */
		index_segment=pointers_to_dists[0]->index;
		target_segment=pointers_to_dists[0]->result;
		hex_join_areas(areas[target_segment],areas[index_segment]);
		areas[index_segment]=NULL;
		for(i=0;i<hx->nareas;++i) {
			if(results[i]==index_segment)
				results[i]=target_segment;
		}
		num_of_removed_areas++;
		pointers_to_dists[0]->distance=2;

		for(i=1;i<num_of_pairs;++i) { /* modify and order indexes */
			pointers_to_dists[i]->modified=0;
			if(pointers_to_dists[i]->index==index_segment) {
				pointers_to_dists[i]->index=target_segment;
				pointers_to_dists[i]->modified=1;
			}

			if(pointers_to_dists[i]->result==index_segment) {
				pointers_to_dists[i]->result=target_segment;
				pointers_to_dists[i]->modified=1;
			}
		}

		for(i=1;i<num_of_pairs;++i) { /* search for areas which are not small*/
			if(pointers_to_dists[i]->modified)  {
				index_segment=pointers_to_dists[i]->index;
				target_segment=pointers_to_dists[i]->result;
				if(index_segment==target_segment) {
					pointers_to_dists[i]->distance=2;
					continue;
				}
				if(!areas[index_segment] || !areas[target_segment]) {
					pointers_to_dists[i]->distance=2;
					continue;
				}
				if(areas[index_segment]->num_of_areas >= pars->minarea) {
					pointers_to_dists[i]->distance=2;
					continue;
				}

				pointers_to_dists[i]->distance=m_find_distance(hx,pars,areas[index_segment],areas[target_segment]);
				for(j=i+1;j<num_of_pairs;++j) { /* search for duplicate entries */
					if(pointers_to_dists[j]->index==pointers_to_dists[i]->index && pointers_to_dists[j]->result==pointers_to_dists[i]->result) {
						pointers_to_dists[j]->distance=2;
						pointers_to_dists[j]->modified=0;
						pointers_to_dists[i]->length_of_edge+=pointers_to_dists[j]->length_of_edge;
					}
				}
			}
		}
		qsort(pointers_to_dists, num_of_pairs, sizeof(DIST*), compare_by_dist);

		if((num_of_removed_areas-previous_num_of_removed_areas)>100)
			fprintf(stderr, "%08d\b\b\b\b\b\b\b\b",num_of_removed_areas);
	x++;
	}
	fprintf(stderr, "\b\b\b\b\b\b\b\b");
	G_message("Number of removed small segments: %d", num_of_removed_areas);

	free(dists);
	free(pointers_to_dists);

	return 0;

}

