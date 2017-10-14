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
#ifdef USE_OPENMP
#include <omp.h>
#endif

struct fifo* hex_get_extended_neighborhood(HEXGRID* hx, struct area** areas, int index)
{
	int i=0,j=0, cur_index1=0,cur_index2=0;
	struct fifo* queue=init_queue(index);

	if(!areas[index])
		return NULL;

	queue=init_queue(index);

	while(hx->hex_neigborhoods[index][i] >-1) {
		cur_index1=hx->hex_neigborhoods[index][i];
		if(!in_queue(queue,cur_index1))
			add_node(queue,cur_index1);
		j=0;
		while(hx->hex_neigborhoods[cur_index1][j] >-1) {
			cur_index2=hx->hex_neigborhoods[cur_index1][j];
			if(!in_queue(queue,cur_index2))
				add_node(queue,cur_index2);
			j++;
		}
		i++;
	}
	pop_node(queue);
	return queue;
}

//################################################

double hex_get_local_heterogeneity2(HEXGRID* hx, LOCAL_PARAMS* p, struct area** areas, int index)
{
	struct fifo* queue;
	int num_of_neighbors;
	int* neighbors;
	int size_of_triangle;
	int* triangle;
	double* distances;
	double mean_distance=0;
	double stddev_distance=0.0;
	int num_of_distances=0;
	double* mean_distances;
	int* mean_distances_num;
	int i,j,k;

	if(!areas[index])
		return -1;

	queue=hex_get_extended_neighborhood(hx,areas,index);
	num_of_neighbors=queue_length(queue);
	neighbors=convert_queue_to_table(queue);

	/* the code works well with triangle */
	size_of_triangle=num_of_neighbors*(num_of_neighbors-1)/2;
	triangle=malloc(size_of_triangle*sizeof(int)); /* crate a vector with indexes of lower trinangle */
	distances=malloc(size_of_triangle*sizeof(double)); /* crate a vector with similarities in lower trinangle */

	/* we create long vector of indexes which pays off to parallelize
	 * in the next step we walk thru the vector and take indexes */

	for(k=0,i=0;i<num_of_neighbors;++i)
		for(j=0;j<i;++j)
			triangle[k++]=i*num_of_neighbors+j;
	if(k!=size_of_triangle)
		G_fatal_error(_("InitSeeds:Wrong_indexing"));

	for(k=0;k<size_of_triangle;++k) { /* reconstruction of matrix */
		int i=triangle[k]/num_of_neighbors;
		int j=triangle[k]%num_of_neighbors;
		double* pair[]={NULL,NULL,NULL};
		if(neighbors[i]>=0)
			pair[0]=hex_use_histogram(hx,neighbors[i]);
		if(neighbors[j]>=0)
			pair[1]=hex_use_histogram(hx,neighbors[j]);
		if(pair[0]==NULL || pair[1]==NULL)
			distances[k]=-1.;
		else
			distances[k]=calculate2(hx,p,pair);
	}
	for(k=0;k<size_of_triangle;++k)
		if(distances[k]>=0) {
			mean_distance+=distances[k];
			stddev_distance+=(distances[k]*distances[k]);
			num_of_distances++;
		}
	mean_distance/=(double)num_of_distances;
	stddev_distance/=(double)num_of_distances;
	stddev_distance=sqrt(stddev_distance-(mean_distance*mean_distance));

	/* now we will try to remove outliers and recalculate mean distance again
	 * we define as outilers those elements which mean distance to other elements is higer than mean_distances + stddev between all elements
	 * similar to inner variance vs outer variance in ANOVA
	 * and define them as outliers */

	mean_distances=malloc(num_of_neighbors*sizeof(double));
	mean_distances_num=calloc(num_of_neighbors,sizeof(int));
	for(i=0;i<num_of_neighbors;++i)
		mean_distances[i]=0.;

	for(k=0;k<size_of_triangle;++k) {
		i=triangle[k]/num_of_neighbors;
		j=triangle[k]%num_of_neighbors;
		if(distances[k]>=0) {
			mean_distances[i]+=distances[k];
			mean_distances[j]+=distances[k];
			mean_distances_num[i]++;
			mean_distances_num[j]++;
		}
	}
	for(i=0;i<num_of_neighbors;++i) {
		if(mean_distances_num[i])
			mean_distances[i]/=mean_distances_num[i];
	}

	int num_of_removed=0;
	for(i=0;i<num_of_neighbors;++i) {
		if(mean_distances[i]>(mean_distance+stddev_distance)) {
			neighbors[i]=-1;
			num_of_removed++;
		}
	}

	mean_distance=0;
	num_of_distances=0;
	for(k=0;k<size_of_triangle;++k) {
		i=triangle[k]/num_of_neighbors;
		j=triangle[k]%num_of_neighbors;
		if(neighbors[i]<0 || neighbors[j]<0 || distances[k]<0)
			continue; /* removed */
		mean_distance+=distances[k];
		stddev_distance+=(distances[k]*distances[k]);
		num_of_distances++;
	}

	mean_distance/=(double)num_of_distances;
	stddev_distance/=(double)num_of_distances;
	stddev_distance=sqrt(stddev_distance-(mean_distance*mean_distance));
	areas[index]->similarity_threshold=mean_distance;//+stddev_distance;
	if(areas[index]->similarity_threshold<0) {
		areas[index]->similarity_threshold=1.;
	}

	free(neighbors);
	free(triangle);
	free(distances);
	free(mean_distances);
	free(mean_distances_num);

	return areas[index]->similarity_threshold;
}

//#####################################################################

double fde(DIST** distances, int cut_point, int num_of_dists)
{
	/* fisher discriminant estimator */
	double mean_distance1=0;
	double mean_distance2=0;
	double var_distance1=0;
	double var_distance2=0;
	int i;

	for(i=0;i<cut_point;++i) {
		mean_distance1+=distances[i]->distance;
		var_distance1+=distances[i]->distance*distances[i]->distance;
	}
	mean_distance1/=((double)cut_point);
	var_distance1/=((double)cut_point);
	var_distance1=var_distance1-(mean_distance1*mean_distance1);

	for(i=cut_point;i<num_of_dists;++i) {
		mean_distance2+=distances[i]->distance;
		var_distance2+=distances[i]->distance*distances[i]->distance;
	}
	mean_distance2/=((double)(num_of_dists-cut_point));
	var_distance2/=((double)(num_of_dists-cut_point));
	var_distance2=var_distance2-(mean_distance2*mean_distance2);

	if(var_distance1==0 && var_distance2==0)
		return -1; /* 100% homogeneous area */

	return ((mean_distance1-mean_distance2)*(mean_distance1-mean_distance2))/
					(var_distance1+var_distance2);
}

int find_fisher_cut(DIST** distances,int num_of_dists, int first)
{

	int i,index=0;
	double estimator, max_estimator=-1;

	for(i=1;i<num_of_dists;++i) {
		estimator=fde(distances,i,num_of_dists);

		if(max_estimator<=estimator) {
			index=i;
			max_estimator=estimator;
		} else {
			if(first)
				return i>2?i-1:2; /* first local maximum */
			}
	}

	return index;
}

double calculate_threshold(DIST** dists, int num_of_neighbors)
{
	int k;
	double mean_distance=0.;
	double stddev_distance=0.;

	for(k=0;k<num_of_neighbors;++k) {
		mean_distance+=dists[k]->distance;
		stddev_distance+=(dists[k]->distance*dists[k]->distance);
	}

	mean_distance/=(double)num_of_neighbors;
	stddev_distance/=(double)num_of_neighbors;
	stddev_distance=sqrt(stddev_distance-(mean_distance*mean_distance));

	return mean_distance+stddev_distance;
}

double hex_get_local_heterogeneity(HEXGRID* hx, LOCAL_PARAMS* p, struct area** areas, int index)
{
	struct fifo* queue;
	int size_of_extended_neighborhood;
	int* extended_neighbors;
	int i;
	DIST* extended_distances;
	DIST** pointers_to_dists;
	double* pair[]={NULL,NULL,NULL};
	int fisher_cut;
	double sequence_position;
	int size_of_direct_neighborhood=0;

	if(!areas[index])
		return -1;

	queue=hex_get_extended_neighborhood(hx,areas,index);
	size_of_extended_neighborhood=queue_length(queue);
	if(size_of_extended_neighborhood<1) {
		areas[index]->similarity_threshold=1;
		free(queue);
		return 3; /* send to the end of queue */
	}
	extended_neighbors=convert_queue_to_table(queue);

	extended_distances=(DIST*)malloc(sizeof(DIST)*size_of_extended_neighborhood);
	pointers_to_dists=(DIST**)malloc(sizeof(DIST*)*size_of_extended_neighborhood);

	pair[0]=hex_use_histogram(hx,index);
	for(i=0;i<size_of_extended_neighborhood;++i) {
		pair[1]=hex_use_histogram(hx,extended_neighbors[i]);
		extended_distances[i].distance=calculate2(hx,p,pair);
		extended_distances[i].index=extended_neighbors[i];
		*(pointers_to_dists+i)=extended_distances+i;
	}

	qsort(pointers_to_dists, size_of_extended_neighborhood, sizeof(DIST*), compare_by_dist);

	fisher_cut=find_fisher_cut(pointers_to_dists,size_of_extended_neighborhood,1); /* first */

	if(fisher_cut==-1) {
		areas[index]->similarity_threshold=0;
	} else {
		areas[index]->similarity_threshold=calculate_threshold(pointers_to_dists,fisher_cut);
		hx->thresholds[index]=calculate_threshold(pointers_to_dists,size_of_extended_neighborhood);
	}

	sequence_position=areas[index]->similarity_threshold;
	size_of_direct_neighborhood=0;
	while(hx->hex_neigborhoods[index][size_of_direct_neighborhood]>-1) size_of_direct_neighborhood++;

	if(size_of_direct_neighborhood<6 || fisher_cut<3)
		sequence_position+=2;

	free(extended_distances);
	free(pointers_to_dists);
	free(extended_neighbors);
	clear_queue(queue);

	return sequence_position;
}

unsigned* hex_find_seeds(HEXGRID* hx, LOCAL_PARAMS* p, struct area** areas, unsigned* num_of_seeds)
{
	int i;
	int j=0;
	int num_of_null_cells=0;
	int* null_flags=calloc(hx->nareas,sizeof(int));
	unsigned* indexes;
	int ncells=hx->nareas;
	DIST* dists=(DIST*)malloc(sizeof(DIST)*ncells);
	DIST** pointers_to_dists=(DIST**)malloc(sizeof(DIST*)*ncells);
	double distance;


	G_message("Finding seeds... (parallel process)");
	Rast_set_d_null_value(hx->thresholds,hx->nareas);

	for(i=0;i<ncells;++i) {
		dists[i].index=i;
		dists[i].distance=1.;
		dists[i].result=(hx->histograms[i]==NULL)?0:1;
		*(pointers_to_dists+i)=dists+i;
	}

#pragma omp parallel for schedule(dynamic,1) private(distance)
	for(i=0;i<ncells;++i) {
		distance=hex_get_local_heterogeneity(hx,p,areas,i);
		if(distance<0)
			null_flags[i]++;
		else {
			dists[i].distance=distance; /* we use distance to estimate sorting position */

			hx->thresholds[i]=areas[i]->similarity_threshold;
			if(areas[i]->similarity_threshold<p->lower_similarity_threshold)
				areas[i]->similarity_threshold=p->lower_similarity_threshold;
			else if(areas[i]->similarity_threshold>p->upper_similarity_threshold)
				areas[i]->similarity_threshold=p->upper_similarity_threshold;
		}
	}
	for(i=0;i<ncells;++i) /* not parrarel */
		num_of_null_cells+=null_flags[i];

	qsort(pointers_to_dists, ncells, sizeof(DIST*), compare_by_dist);
	indexes=malloc((ncells-num_of_null_cells)*sizeof(unsigned));

	for(i=0;i<ncells;++i) {
		if(pointers_to_dists[i]->result) {
			indexes[j++]=pointers_to_dists[i]->index;
		}
	}

	/* free memory */
	free(pointers_to_dists);
	free(dists);
	free(null_flags);
	/* return */
	*num_of_seeds=j;

	return indexes;
}


double* create_thresholds_map(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas)
{
	/* convert to regular grid */
	int nrows=d->cell_hd.rows;
	int ncols=d->cell_hd.cols;
	int ncells=nrows*ncols;
	int i,k;
	double* results=malloc(ncells*sizeof(double));
	int* histogram_ids=NULL;

	Rast_set_d_null_value(results,ncells);

	for(i=0;i<hx->nareas;++i) {
		histogram_ids=hx->histogram_ids[i];
		if(histogram_ids==NULL) /* can happen */
			continue;

		for(k=0;k<hx->size_of_supermotifel;++k) {
			if(histogram_ids[k]>-1)
				results[histogram_ids[k]]=hx->thresholds[i];
		}
	}

	return results;
}
