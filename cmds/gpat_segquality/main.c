/****************************************************************************
 *
 * PROGRAM:    gpat_segquality - part of GeoPAT 2
 * AUTHOR(S):   Jacek Niesterowicz, Pawel Netzel
 * PURPOSE:    program to calculate quality of segmentation;
 *        segmentation code by Jaroslaw Jasiewicz, Jacek Niesterowicz, 
 *        Tomasz Stepinski from GRASS GeoPAT 
 * COPYRIGHT:    (C) Jacek Niesterowicz, Pawel Netzel, Space Informatics Lab,
 *        University of Cincinnati
 *              http://sil.uc.edu
 *
 *        This program is free software under 
 *        the GNU General Public License (>=v3). 
 *        https://www.gnu.org/licenses/gpl-3.0.en.html
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
    struct area** areas=NULL;
    HEXGRID* hexgrid;
    DATAINFO* datainfo;
    LOCAL_PARAMS* parameters;
    parameters=malloc(sizeof(LOCAL_PARAMS));
    datainfo=malloc(sizeof(DATAINFO));
    char *list_dist;
    double* heterogeneity_map;
    double* isolation_map;
    int* segment_map;
    char *list;
    int size_val = 1;
    EZGDAL_LAYER *input_layer;

    struct arg_str  *inp   = arg_str1("i","input","<file_name>","name of input file (GRID)");
    struct arg_str  *seg   = arg_str1("s","segments","<file name>","name of input segmentation map (TIFF)");
    struct arg_str  *het   = arg_str0("g","inhomogeneity","<file_name>","name of output file with segment inhomogeneity (TIFF)");
    struct arg_str  *iso   = arg_str0("o","isolation","<file_name>","name of output file with segment isolation (TIFF)");
    struct arg_str  *mes   = arg_str0("m","measure","<measure_name>","similarity measure (use -l to list all measures; defult: jsd)");
    struct arg_lit  *mesl  = arg_lit0("l","list_measures","list all measures");
    struct arg_int  *maxhist       = arg_int0(NULL,"maxhist","<n>","create similarity/distance matrix for maxhist histograms; leave 0 to use all (default: 200)");
    struct arg_lit  *flag_complete = arg_lit0("c","complete","use complete linkage (default is average)");
    struct arg_lit  *flag_quad     = arg_lit0("q","quad","quad mode (rook topology)");
    struct arg_lit  *flag_noweight = arg_lit0("w","no_weight","switch off edge-based weighting in isolation");
    struct arg_int  *th            = arg_int0("t",NULL,"<n>","number of threads (default: 1)");
    struct arg_lit  *help          = arg_lit0("h","help","print help and exit");
    struct arg_end  *end           = arg_end(20);

    void* argtable[] = {inp,seg,het,iso,mes,mesl,
                        maxhist,flag_complete,flag_quad,
                        flag_noweight,th,help,end};

    int nerrors = arg_parse(argc,argv,argtable);

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

    if(het->count==0 && iso->count==0){
      printf("\nAt least one of the outputs (inhomogeneity, isolation) must be provided.\n\n");
      exit(0);
    }

    if(th->count > 0) 
      omp_set_num_threads(th->ival[0]);
    else
      omp_set_num_threads(1);

    if(!ezgdal_file_exists(inp->sval[0])) {
      printf("\nInput file [%s] does not exist\n\n",inp->sval[0]);
      exit(0);
    }

    if(!ezgdal_file_exists(seg->sval[0])) {
      printf("\nInput map [%s] does not exist\n\n",seg->sval[0]);
      exit(0);
    }
    
	/* parameters preparation */

	if(maxhist->count>0)
      parameters->sampling_threshold=maxhist->ival[0];
    else
      parameters->sampling_threshold=200;
    if(parameters->sampling_threshold<0) {
      printf("\nSampling threshold cannot be negative\n\n");
      exit(0);
    }

    parameters->reduction=0;

    parameters->quad_mode=(flag_quad->count>0);
    parameters->complete_linkage=(flag_complete->count>0);
    parameters->no_weight=(flag_noweight->count>0);
    parameters->null_threshold = 0.5;

    if(mes->count > 0) {
      parameters->calculate = get_distance((char *)(mes->sval[0]));
      /* measure not found */
      if(parameters->calculate==NULL) {
        printf("\nWrong distance measure: %s\n\n",mes->sval[0]);
        printf("list of distances:\n");
        list_dist = list_all_distances();
        printf("\n%s\n",list_dist);
        free(list_dist);
        exit(0);
      }
    } else 
      parameters->calculate = get_distance("jsd");

	/* open segmentation map */
    input_layer = (EZGDAL_LAYER *)malloc(sizeof(EZGDAL_LAYER));
    input_layer = ezgdal_open_layer((char *)(seg->sval[0]));
    if(input_layer==NULL) {
      printf("\nCannot open file: '%s'\n\n", seg->sval[0]);
    }

    printf("Preparing..."); fflush(stdout);

	/* open and read grid */
    datainfo=malloc(sizeof(DATAINFO));
    init_grid_datainfo(datainfo,(char *)(inp->sval[0]),"apud");
    read_signatures_to_memory(datainfo);
    parameters->parameters = init_measure_parameters(datainfo->size_of_histogram,0); /* use distance, not similarity */
    
	/* check extents */
	if(datainfo->cell_hd.rows!=input_layer->rows || datainfo->cell_hd.cols!=input_layer->cols) {
      printf("\nSize of input grid and map do not match.\n\n");
      exit(0);
	}

	/* rebuild segmentation topology */
    hexgrid = hex_build_topology(datainfo,parameters);
	segment_map = create_segment_map(hexgrid, parameters, input_layer);
	areas=build_areas(datainfo,hexgrid,parameters,segment_map);

    printf("OK\n"); fflush(stdout);

	/* calculate diagnostics */
    if(het->count > 0) {
        heterogeneity_map = calculate_heterogeneity(datainfo, hexgrid, parameters, areas);
        write_raster_dbl(datainfo->dh, (void*)heterogeneity_map, (char *)(het->sval[0]), size_val);
        free(heterogeneity_map);
    }
    if(iso->count > 0) {
		build_neighbors_queue(hexgrid,parameters,areas,segment_map);
	    isolation_map = calculate_isolation(datainfo, hexgrid, parameters, areas);
        write_raster_dbl(datainfo->dh, (void*)isolation_map, (char *)(iso->sval[0]), size_val);
        free(isolation_map);
    }

	/* free */
	ezgdal_close_layer(input_layer);
	free(segment_map);
	remove_areas(areas,hexgrid->nareas);
    hex_remove_hexgrid(hexgrid);
    free(datainfo);
    free(parameters);
    free(hexgrid);

    exit(0);
}

