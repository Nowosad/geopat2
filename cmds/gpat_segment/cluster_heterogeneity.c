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
#ifdef USE_OPENMP
#include <omp.h>
#endif

/*========================================= homogenity ======================================*/
double hex_calculate_linkage(HEXGRID* hx, LOCAL_PARAMS* p, struct area* a, struct area* b)
{
	double* pair[]={NULL,NULL,NULL};
	int index0;
	int index1;
	int i,j,k;
	int *a_samples=NULL, *b_samples=NULL;
	int num_of_a_samples, num_of_b_samples;
	int size_of_matrix;
	int* matrix_indexes;
	double* distances;
	int* null_flags;
	int num_of_null_pairs=0;
	double linkage=0;

	if(a->hex_ids==NULL)
		a->hex_ids=init_queue(a->id);

	if(b->hex_ids==NULL)
		b->hex_ids=init_queue(b->id);

	if(a->num_of_areas==1 && b->num_of_areas==1) {
		index0=get_value_at_position(a->hex_ids,0);
		index1=get_value_at_position(b->hex_ids,0);
		pair[0]=hex_use_histogram(hx,index0);
		pair[1]=hex_use_histogram(hx,index1);

		return calculate2(hx,p,pair);
	}

	/* sampling */
	if(a->num_of_areas<=(int)p->sampling_threshold*1.1||p->sampling_threshold==0)
		num_of_a_samples=a->num_of_areas;
	else
		num_of_a_samples=p->sampling_threshold;
	a_samples=sample_histogram_ids(a->hex_ids,num_of_a_samples);

	if(b->num_of_areas<=(int)p->sampling_threshold*1.1||p->sampling_threshold==0)
		num_of_b_samples=b->num_of_areas;
	else
		num_of_b_samples=p->sampling_threshold;
	b_samples=sample_histogram_ids(b->hex_ids,num_of_b_samples);
	/* end of sampling */

	size_of_matrix=num_of_a_samples*num_of_b_samples;
	/* matrix is asimmetrical a are arannged as rows, b are arranged as cols
	 * index is calculated as a * num_of_b_samples + b
	 * to retrive col and row from index use
	 * a/num_of_b_samples for rows and
	 * a%num_of_b_samples for cols and */

	/* this approach is to make parrarel code efficient */
	matrix_indexes=malloc(size_of_matrix*sizeof(int));
	for(i=0;i<size_of_matrix;++i)
		matrix_indexes[i]=i;

	distances=malloc(size_of_matrix*sizeof(double));
	null_flags=calloc(size_of_matrix,sizeof(double));

#pragma omp parallel for schedule(dynamic,1) private(i,j,pair)
	for(k=0;k<size_of_matrix;++k) { /* reconstruction of matrix */
		i=matrix_indexes[k]/num_of_b_samples;
		j=matrix_indexes[k]%num_of_b_samples;
		pair[0]=NULL;
		pair[1]=NULL;
		pair[0]=hex_use_histogram(hx,a_samples[i]);
		pair[1]=hex_use_histogram(hx,b_samples[j]);
		if(pair[0]==NULL || pair[1]==NULL) {
			distances[k]=0.;
			null_flags[k]=0;
		}
		else
			distances[k]=calculate2(hx,p,pair);
	}

	free(a_samples);
	free(b_samples);
	free(matrix_indexes);

	if(p->complete_linkage) {
		for(k=0;k<size_of_matrix;++k) {
			if(distances[k]>linkage)
				linkage=distances[k];
		}
	}
	else {
		for(k=0;k<size_of_matrix;++k) {
			linkage+=distances[k];
			num_of_null_pairs+=null_flags[k];
		}
		if((size_of_matrix-num_of_null_pairs)==0) {
			free(distances);
			free(null_flags);
			return 1.;
		}
		linkage/=(double)(size_of_matrix-num_of_null_pairs);
	}
	free(distances);
	free(null_flags);

	return linkage;

	/* return common homogeneity */
//	int num_a_dist=a->num_of_areas*(a->num_of_areas-1)*0.5;
//	int num_b_dist=b->num_of_areas*(b->num_of_areas-1)*0.5;
//	double value=(a->homogenity*num_a_dist+b->homogenity*num_b_dist+linkage)/(num_a_dist+num_b_dist+size_of_matrix-num_of_null_pairs);
//	printf("%d %d, %f \n",a->id,b->id,value);
//	return value;

}

int add_areas(struct area* a, struct area* b)
{
	return a->num_of_areas+b->num_of_areas;
}

int hex_join_areas(struct area* a, struct area* b)
{
	join_queues(a->hex_ids,b->hex_ids);
	a->num_of_areas=add_areas(a,b);
	free(b->hex_ids);
	b->hex_ids=NULL;

	remove_area(&b);
	b=NULL;

	return 0;
}
