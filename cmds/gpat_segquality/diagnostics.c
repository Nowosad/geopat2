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

/* =========================== HETEROGENEITY ================================ */

int calculate_segment_heterogeneity(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* pars, struct area* a)
{
	int i,j;
    double mean=0.;
    double* pair[]={NULL,NULL,NULL};
    int num_of_pairs=0;
	distance_func* calculate=pars->calculate;
	int dims[1];
	dims[0]=pars->parameters->size_of_histogram;

    if(a->num_of_samples<2) { /* nothing to calculate */
        a->heterogeneity=0;
        return 0;
    }
//printf("%d\n",a->id);
    for(i=0;i<a->num_of_samples;++i) {
        for(j=(i+1);j<a->num_of_samples;++j) {
            pair[0]=hex_use_histogram(hx,a->sample_ids[i]);
            pair[1]=hex_use_histogram(hx,a->sample_ids[j]);
            mean+=calculate(pair,2,pars->parameters->size_of_histogram,1,dims);
            num_of_pairs++;
        }
    }
    a->heterogeneity=mean/(double)num_of_pairs;

    return 0;
}

/* main function */
double* calculate_heterogeneity(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas)
{
    int i,j,k,c=0;
    int num_of_areas=0;
    int* histogram_ids=NULL;
    int* ids=NULL;
    int num_of_histograms;
    double* heterogeneity_map;

    printf("Calculating inhomogeneity...\n");

    heterogeneity_map = malloc(d->cell_hd.rows * d->cell_hd.cols * sizeof(double));
	Rast_set_d_null_value(heterogeneity_map,d->cell_hd.rows * d->cell_hd.cols);

    for(i=0;i<hx->nareas;++i) {
        if(areas[i])
            num_of_areas++;
    }
    for(i=0;i<hx->nareas;++i) {
        if(areas[i]==NULL)
        	continue;
        ezgdal_show_progress(stdout,c++,num_of_areas);
        calculate_segment_heterogeneity(d, hx, pars, areas[i]);

        /* write to map */
	    ids=areas[i]->hex_ids;
	    num_of_histograms=areas[i]->num_of_areas;
	    for (j=0;j<num_of_histograms;++j) {
	        histogram_ids=hx->histogram_ids[ids[j]];
	        for(k=0;k<hx->size_of_supermotifel;++k) {
	            if(histogram_ids[k]>-1) {
	                heterogeneity_map[histogram_ids[k]]=areas[i]->heterogeneity;
	            }
	        }
	    }
    }
    ezgdal_show_progress(stdout,100,100);

    return heterogeneity_map;
}

/* ============================== ISOLATION ================================= */

double linkage(HEXGRID* hx, LOCAL_PARAMS* p, struct area* a, struct area* b)
{
    double* pair[]={NULL,NULL,NULL};
    int index0;
    int index1;
    int i,j,k;
    int size_of_matrix;
    int* matrix_indexes;
    double* distances;
    int* null_flags;
    int num_of_null_pairs=0;
    double linkage=0;
	distance_func* calculate=p->calculate;
	int dims[1];
	dims[0]=p->parameters->size_of_histogram;
	struct fifo* hex_ids;

    if(a->hex_ids==NULL){
        hex_ids=init_queue(a->id);
        a->hex_ids=convert_queue_to_table(hex_ids);
		clear_queue(hex_ids);
        free(hex_ids);
    }
    if(b->hex_ids==NULL){
        hex_ids=init_queue(a->id);
        b->hex_ids=convert_queue_to_table(hex_ids);
		clear_queue(hex_ids);
        free(hex_ids);
	}

    if(a->num_of_areas==1 && b->num_of_areas==1) {
        index0=a->hex_ids[0];
        index1=b->hex_ids[0];
        pair[0]=hex_use_histogram(hx,index0);
        pair[1]=hex_use_histogram(hx,index1);

        return calculate(pair,2,p->parameters->size_of_histogram,1,dims);
    }

    size_of_matrix=a->num_of_samples * b->num_of_samples;
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
        i=matrix_indexes[k]/b->num_of_samples;
        j=matrix_indexes[k]%b->num_of_samples;
        pair[0]=NULL;
        pair[1]=NULL;
        pair[0]=hex_use_histogram(hx,a->sample_ids[i]);
        pair[1]=hex_use_histogram(hx,b->sample_ids[j]);
        if(pair[0]==NULL || pair[1]==NULL) {
            distances[k]=0.;
            null_flags[k]=0;
        }
        else
            distances[k]=calculate(pair,2,p->parameters->size_of_histogram,1,dims);
    }

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
}

int calculate_segment_isolation(HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas, struct area* a)
{
    int i;
    int sum_length_of_edges=0;
    double mean_dist=0.;
    int target;

    for(i=1;i<a->num_of_neighbors;++i) { /* first element contains id of the segment shall we check? */
        target=a->neighbors[i];
        double lk=linkage(hx,pars,a,areas[target]);
        if(pars->no_weight){
    		mean_dist+=lk;
    		sum_length_of_edges+=1;
    	}
    	else{
	        mean_dist+=(lk*a->edges[i]);
    	    sum_length_of_edges+=a->edges[i];
    	}
    }
    a->isolation=mean_dist/(double)sum_length_of_edges;

    return 0;
}

/* main function */
double* calculate_isolation(DATAINFO* d, HEXGRID* hx, LOCAL_PARAMS* pars, struct area** areas)
{
    int i,j,k,c=0;
    int num_of_areas=0;
    int* histogram_ids=NULL;
    int* ids=NULL;
    int num_of_histograms;
    double* isolation_map;

    printf("Calculating isolation...\n");

    isolation_map = malloc(d->cell_hd.rows * d->cell_hd.cols * sizeof(double));
	Rast_set_d_null_value(isolation_map, d->cell_hd.rows * d->cell_hd.cols);

    for(i=0;i<hx->nareas;++i) {
        if(areas[i])
            num_of_areas++;
    }
    for(i=0;i<hx->nareas;++i) {
        if(areas[i]==NULL)
        	continue;
        ezgdal_show_progress(stdout,c++,num_of_areas);
        calculate_segment_isolation(hx, pars, areas, areas[i]);

        /* write to map */
        ids=areas[i]->hex_ids;
        num_of_histograms=areas[i]->num_of_areas;
	    for (j=0;j<num_of_histograms;++j) {
	        histogram_ids=hx->histogram_ids[ids[j]];
	        for(k=0;k<hx->size_of_supermotifel;++k) {
	            if(histogram_ids[k]>-1) {
	                isolation_map[histogram_ids[k]]=areas[i]->isolation;
	            }
	        }
	    }
    }
    ezgdal_show_progress(stdout,100,100);

    return isolation_map;
}

