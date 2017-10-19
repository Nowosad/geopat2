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
double find_distance(HEXGRID* hx, LOCAL_PARAMS* pars, struct area* index, struct area* target)
{
	double tmp_dist=hex_calculate_linkage(hx,pars,index,target);
	double min_similarity=MAX(pars->lower_similarity_threshold,MIN(index->similarity_threshold,target->similarity_threshold));

	return tmp_dist-min_similarity;
}

int hex_hierarhical(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, int *results)
{
	int i,j,num_of_pairs=0;
	int index_segment,target_segment,tmp_index;
	int num_of_removed_areas=0;
	double similarity_threshold;
	DIST* dists;
	DIST** pointers_to_dists;
	int k=0; /* index of pairs */

	G_message("Removing segments by joining...");
	for(i=0;i<hx->nareas;++i) {
		if(areas[i]!=NULL) {
			if(areas[i]->neighbors!=NULL) {
				clear_queue(areas[i]->neighbors); /* clear all areas which contain queue */
				areas[i]->neighbors=NULL;
				areas[i]->num_of_neighbors=0;
			}
			areas[i]->neighbors=init_queue(i); /* first element of queue is area, next elements will be added after it! */
			areas[i]->done=0;
		}
	}
	fprintf(stderr, "...Buliding graph...\n");
	hex_build_neighbors_queue(hx,pars,areas,results,1); /* remove  dual pairs */

	for(i=0;i<hx->nareas;++i)
		if(areas[i]!=NULL && areas[i]->neighbors->length>1) /* first element is an index, so one element means segments has no joins */
			num_of_pairs+=(areas[i]->neighbors->length-1);

	dists=(DIST*)malloc(sizeof(DIST)*num_of_pairs);
	pointers_to_dists=(DIST**)malloc(sizeof(DIST*)*num_of_pairs);

	/* convert list into table, maybe not very clear */
	for(i=0;i<hx->nareas;++i) {
		if(areas[i]!=NULL && areas[i]->neighbors->length>0) { /* as long as there are some neighbors */
			index_segment=pop_node(areas[i]->neighbors);
			while(areas[index_segment]->neighbors->first) {
				target_segment=pop_node(areas[index_segment]->neighbors); /* we remove all targets */
				dists[k].index=index_segment;
				dists[k].result=target_segment; /* use results as targets */
				dists[k].distance=find_distance(hx,pars,areas[index_segment],areas[target_segment]);
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
	fprintf(stderr, " %d pairs Done \n", num_of_pairs);

	while (pointers_to_dists[0]->distance<=0) {
		/* areas are sorted ASC by distance all others have distance grater than threshold so we can break for loop */
		index_segment=pointers_to_dists[0]->index;
		target_segment=pointers_to_dists[0]->result;
		similarity_threshold=MIN(areas[index_segment]->similarity_threshold,areas[target_segment]->similarity_threshold);
		hex_join_areas(areas[index_segment],areas[target_segment]);
		areas[index_segment]->similarity_threshold=similarity_threshold;
		areas[target_segment]=NULL;
		for(i=0;i<hx->nareas;++i)
			if(results[i]==target_segment)
				results[i]=index_segment;
		num_of_removed_areas++;
		pointers_to_dists[0]->distance=1;

		for(i=1;i<num_of_pairs;++i) { /* modify and order indexes */
			pointers_to_dists[i]->modified=0;
			if(pointers_to_dists[i]->index==target_segment) {
				pointers_to_dists[i]->index=index_segment;
				pointers_to_dists[i]->modified=1;
			}

			if(pointers_to_dists[i]->result==target_segment) {
				pointers_to_dists[i]->result=index_segment;
				pointers_to_dists[i]->modified=1;
			}

			if(pointers_to_dists[i]->modified && pointers_to_dists[i]->result < pointers_to_dists[i]->index) { /* swap to have ascending order */
				tmp_index=pointers_to_dists[i]->result;
				pointers_to_dists[i]->result=pointers_to_dists[i]->index;
				pointers_to_dists[i]->index=tmp_index;
			}
		}

		for(i=1;i<num_of_pairs;++i) { /* search for duplicate entries */
			if(pointers_to_dists[i]->modified)  {
				index_segment=pointers_to_dists[i]->index;
				target_segment=pointers_to_dists[i]->result;
				if(index_segment==target_segment) {
					pointers_to_dists[i]->distance=2;
					continue;
				}
				if(!areas[index_segment] || !areas[target_segment]) {
					pointers_to_dists[i]->distance=1;
					continue;
				}
				pointers_to_dists[i]->distance=find_distance(hx,pars,areas[index_segment],areas[target_segment]);
				for(j=i+1;j<num_of_pairs;++j) {
					if(pointers_to_dists[j]->index==pointers_to_dists[i]->index && pointers_to_dists[j]->result==pointers_to_dists[i]->result) {
						pointers_to_dists[j]->distance=1;
						pointers_to_dists[j]->modified=0;
					}
				}
			}
		}
		qsort(pointers_to_dists, num_of_pairs, sizeof(DIST*), compare_by_dist);
		fprintf(stderr, "%08d\b\b\b\b\b\b\b\b",num_of_removed_areas);
	}
	fprintf(stderr, "\b\b\b\b\b\b\b\b");

	free(dists);
	free(pointers_to_dists);
	G_message("Number of removed segments: %d", num_of_removed_areas);

	return 0;

}
