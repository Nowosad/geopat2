/****************************************************************************
 *
 * PROGRAM:	gpat_segment - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel, Jakub Nowosad
 * PURPOSE:	program for a grid of motifels segmentation;
 *		segmentation code by Jaroslaw Jasiewicz, Jacek Niesterowicz, 
 *		Tomasz Stepinski from GRASS GeoPAT 
 * COPYRIGHT:	(C) Pawel Netzel, Space Informatics Lab,
 *		University of Cincinnati
 *              http://sil.uc.edu
 *
 *		This program is free software under 
 *		the GNU General Public License (>=v3). 
 *		https://www.gnu.org/licenses/gpl-3.0.en.html
 *
 *****************************************************************************/
#define MAIN
#include <omp.h>
#include "local_proto.h"

#include "../../lib/ezGDAL/ezgdal.h"

#include "../../lib/argtable/argtable3.h"
#include "../../lib/measures/measures.h"
#include "../../lib/tools/libtools.h"

void usage(char *progname, void *argtable) {
      printf("\nUsage:\n\t%s", progname);
      arg_print_syntax(stdout,argtable,"\n");
      printf("\n");
      arg_print_glossary_gnu(stdout,argtable);
      printf("\n");
      exit(0);
}

int main(int argc, char *argv[])
{
//    Cell_head window;
    CELL* results=NULL;
    struct area** areas=NULL;
    HEXGRID* hexgrid;
    DATAINFO** datainfo;
    LOCAL_PARAMS* parameters;
    parameters=malloc(sizeof(LOCAL_PARAMS));
    datainfo=malloc(sizeof(DATAINFO));
    char *list_dist;
    unsigned num_of_seeds;
    unsigned* seeds;
    int* segment_map;
    char *list;
    int size_val = 1;
    char *weights_val = NULL;

    struct arg_str  *inp   = arg_strn("i","input","<file_name>",1,9999,"name of input files (GRID)");
    struct arg_str  *out   = arg_str1("o","output","<file_name>","name of output file with segments (TIFF)");
    struct arg_str  *shp   = arg_str0("v","vector","<file_name>","name of output vector file with segments (GPKG)");
    struct arg_int  *size  = arg_int0(NULL,"size","<n>","output resolution modifier (default: 1)");
    struct arg_str  *mes   = arg_str0("m","measure","<measure_name>","similarity measure (use -l to list all measures; default: jsd)");
    struct arg_lit  *mesl  = arg_lit0("l","list_measures","list all measures");
    struct arg_dbl  *lower_threshold  = arg_dbl0(NULL,"lthreshold","<double>","minimum distance threshold to build areas (default: 0.1)");
    struct arg_dbl  *upper_threshold  = arg_dbl0(NULL,"uthreshold","<double>","maximum distance threshold to build areas (default: 0.3");
    struct arg_str  *weights    = arg_str0("w","weights","<integer,integer>","multilayer only: weights for the multilayer mode");
    struct arg_dbl  *swap       = arg_dbl0(NULL,"swap","<double>","improve segmentation by swapping unmatched areas. -1 to skip (default: 0.001)");
    struct arg_int  *minarea    = arg_int0(NULL,"minarea","<n>","minimum number of motifels in individual segment (default: 0)");
    struct arg_int  *maxhist    = arg_int0(NULL,"maxhist","<n>","create similarity/distance matrix for maxhist histograms; leave 0 to use all (default: 200)");
    struct arg_lit  *flag_complete          = arg_lit0("c","complete","use complete linkage (default is average)");
//  struct arg_lit  *flag_threshold         = arg_lit0("d","th_map","calculate threshold layer and exit (all params are ignored)");
    struct arg_lit  *flag_skip_growing      = arg_lit0("g","no_growing","skip growing phase");
    struct arg_lit  *flag_skip_hierarchical = arg_lit0("r","no_hierarchical","skip hierarchical phase");
//    struct arg_lit  *flag_all               = arg_lit0("a","all_layers","multilayer only: compare a threshold against all layers instead of an average");
    struct arg_lit  *flag_quad              = arg_lit0("q","quad","quad mode (rook topology)");
    struct arg_int  *th    = arg_int0("t",NULL,"<n>","number of threads (default: 1)");
    struct arg_lit  *help  = arg_lit0("h","help","print help and exit");
    struct arg_end  *end   = arg_end(20);

    void* argtable[] = {inp,out,shp,size,mes,mesl,
                        lower_threshold,upper_threshold,weights,
                        swap,minarea,maxhist,
                        flag_complete,/*flag_threshold,*/flag_skip_growing,
                        flag_skip_hierarchical,/*flag_all,*/flag_quad,
                        th,help,end};

    int nerrors = arg_parse(argc,argv,argtable);
    int num_of_layers = inp->count;

    if (help->count > 0) 
      usage(argv[0],argtable);

    /* list all measures */
    if(mesl->count > 0) {
      list = list_all_distances();
      printf("\nList of measures:\n\n%s\n",list);
      free(list);
      exit(0);
    }
    
    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0) {
      /* Display the error details contained in the arg_end struct.*/
      arg_print_errors(stdout,end,argv[0]);
      usage(argv[0],argtable);
    }

    if(th->count > 0) 
      omp_set_num_threads(th->ival[0]);
    else
      omp_set_num_threads(1);
  
    int i;
    for(i=0; i<inp->count; i++){
      if(!ezgdal_file_exists(inp->sval[i])) {
        printf("\nInput file [%s] does not exist\n\n",inp->sval[i]);
      usage(argv[0],argtable);
      }
    }

    if(flag_skip_growing->count>0 && flag_skip_hierarchical->count>0) {
      printf("\nOnly one flag of -r or -g can be used\n\n");
      exit(0);
    }

    if(lower_threshold->count>0)
      parameters->lower_similarity_threshold=lower_threshold->dval[0];
    else
      parameters->lower_similarity_threshold=0.1;
    if(parameters->lower_similarity_threshold<0.0 || parameters->lower_similarity_threshold>1.0) {
      printf("\nLower distance threshold must be between 0 and 1\n\n");
      exit(0);
    }

    if(size->count>0) {
      size_val = size->ival[0];
      if(size_val<=0) size_val = 1;
    }

    if(upper_threshold->count>0)
      parameters->upper_similarity_threshold=upper_threshold->dval[0];
    else
      parameters->upper_similarity_threshold=0.3;
    if(parameters->upper_similarity_threshold<0.0 || parameters->upper_similarity_threshold>1.0) {
      printf("\nUpper distance threshold must be between 0 and 1\n\n");
      exit(0);
    }
    if(parameters->upper_similarity_threshold < parameters->lower_similarity_threshold) {
      printf("\nUpper distance threshold cannot be smaller than lower threshold\n\n");
      exit(0);
    }
    
    if(weights->count==0){
      parameters->all_layers = 1;
    } else if (weights->count>0){
      parameters->all_layers = 0;
      weights_val = (char *)(weights->sval[0]);
    }
    
    if(swap->count>0)
      parameters->swap_threshold=swap->dval[0];
    else
      parameters->swap_threshold=0.001;
    if(parameters->swap_threshold>1.0) {
      printf("\nSwap must be between 0 and 1 or negative to skip\n\n");
      exit(0);
    }

    if(maxhist->count>0)
      parameters->sampling_threshold=maxhist->ival[0];
    else
      parameters->sampling_threshold=200;
    if(parameters->sampling_threshold<0) {
      printf("\nSampling threshold cannot be negative\n\n");
      exit(0);
    }

    if(minarea->count>0)
      parameters->minarea=minarea->ival[0];
    else
      parameters->minarea=0;
    if(parameters->minarea<0) {
      printf("\nArea must be non-negative\n\n");
      exit(0);
    }

    parameters->reduction=0;

    parameters->quad_mode=(flag_quad->count>0);
    parameters->complete_linkage=(flag_complete->count>0);
    parameters->null_threshold = 0.5;
//    parameters->all_layers=(flag_all->count>0);
    
//    if(parameters->all_layers && weights_val)
//        G_warning("Ignore weigths in the <all layers> mode");

    if(mes->count > 0) {
      parameters->calculate = get_distance((char *)(mes->sval[0]));
      /* measure not found */
      if(parameters->calculate==NULL) {
        printf("\nWrong distance measure: %s\n\n",mes->sval[0]);
        printf("List of distances:\n");
        list_dist = list_all_distances();
        printf("\n%s\n",list_dist);
        free(list_dist);
        exit(0);
      }
    } else 
      parameters->calculate = get_distance("jsd");

	  datainfo = malloc(num_of_layers*sizeof(DATAINFO*));

	  
	  for(i=0; i<num_of_layers; ++i) {
	    datainfo[i] = malloc(num_of_layers*sizeof(DATAINFO*));
      init_grid_datainfo(datainfo[i],(char *)(inp->sval[i]),(char *)(out->sval[0]));

      read_signatures_to_memory(datainfo[i]);
    }

    hexgrid = hex_build_topology(datainfo,parameters,num_of_layers,weights_val);
    areas = hex_build_areas(datainfo,hexgrid,parameters);
    results = hex_init_results(hexgrid);
    parameters->parameters = init_measure_parameters(datainfo[0]->size_of_histogram,0); /* we will use distance instead of similarity */

   /* seeding starts here */
    seeds = hex_find_seeds(hexgrid,parameters,areas,&num_of_seeds);

  /* thresholds only */
/*
	if(flag_threshold->answer!=0) {
		double* thresholds=create_thresholds_map(datainfo[0],hexgrid,parameters,areas);
		Rast_get_window(&window);
		G_message(_("Change window to write results"));
		Rast_set_window(&(datainfo[0]->cell_hd));
		char map_name[100]="\0";
		sprintf(map_name,"%s_threshold",opt_output_layer->answer);
		write_map(datainfo[0],parameters,(void*)thresholds,map_name,DCELL_TYPE,2);
		G_message(_("Restore original region definition"));
		Rast_set_window(&window);
		free(seeds);
		hex_remove_hexgrid(hexgrid);
		free(results);
		free(datainfo);
		free(parameters);
		free(hexgrid);
		exit(EXIT_SUCCESS);
	}
*/

    /* segmentation starts here */
    if(flag_skip_growing->count==0)
        hex_region_growing(hexgrid,parameters,areas,results,seeds,num_of_seeds);

    if(flag_skip_hierarchical->count==0)
        hex_hierarchical(hexgrid,parameters,areas,results);

    if(parameters->minarea>0)
        hex_minarea(hexgrid,parameters,areas,results);

    if(parameters->swap_threshold>=0) {
         swap_areas(hexgrid,parameters,areas,results);
         if(parameters->minarea>0)
             hex_minarea(hexgrid,parameters,areas,results);
    }

    hex_reclass(hexgrid,areas);
    segment_map=hex_create_segment_map(datainfo[0],hexgrid,parameters,areas);

    write_raster(datainfo[0]->dh, (void*)segment_map, (char *)(out->sval[0]), size_val);
    if(shp->count>0)
	convert_to_vector((char *)(out->sval[0]),(char *)(shp->sval[0]));
	
    hex_remove_hexgrid(hexgrid);
    free(segment_map);
    free(results);
    free(datainfo);
    free(parameters); /* TODO: ADD RELEASE FUNCTION */;
    free(hexgrid);


    exit(0);
}

