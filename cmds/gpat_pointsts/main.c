/****************************************************************************
 *
 * PROGRAM:	gpat_pointsts - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	program for a time serie of given point or a set of points;
 * COPYRIGHT:	(C) Pawel Netzel, Space Informatics Lab,
 *		University of Cincinnati
 *              http://sil.uc.edu
 *
 *		This program is free software under 
 *		the GNU General Public License (>=v3). 
 *		https://www.gnu.org/licenses/gpl-3.0.en.html
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include <sml.h>
#include <ezgdal.h>

#include "../../lib/argtable/argtable3.h"
#include "../../lib/tools/libtools.h"

void usage(char *progname, void *argtable) {
      printf("\nUsage:\n\t%s", progname);
      arg_print_syntax(stdout,argtable,"\n");
      printf("\n");
      arg_print_glossary_gnu(stdout,argtable);
      printf("\n");
      exit(0);
}

int main(int argc, char **argv) {

    char *description;
    SML_DATA_HEADER *dh_f;
    void *xybuf;
    FILE *f, *fxy;
    char desc_text[MAX_DESC_LEN];
    double coord_x,coord_y;

    struct arg_str  *inp  = arg_str1("i","input","<file_name>","name of input file (GRID)");
    struct arg_str  *out  = arg_str1("o","output","<file_name>","name of output file (TXT)");
    struct arg_dbl  *x    = arg_dbl0("x",NULL,"<double>","x coord");
    struct arg_dbl  *y    = arg_dbl0("y",NULL,"<double>","y coord");
    struct arg_str  *desc = arg_str0("d","description","<string>","Description of the location");
    struct arg_str  *xy   = arg_str0(NULL,"xy_file","<file_name>","name of file with coordinates (TXT)");
    struct arg_lit  *app  = arg_lit0("a","append","append results to output file");
    struct arg_lit  *help = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end  = arg_end(20);
    void* argtable[] = {inp,out,x,y,desc,xy,app,help,end};

    int nerrors = arg_parse(argc,argv,argtable);

    description = "location";

    if (help->count > 0) 
      usage(argv[0],argtable);

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0) {
      /* Display the error details contained in the arg_end struct.*/
      arg_print_errors(stdout,end,argv[0]);
      usage(argv[0],argtable);
    }

    if(!((x->count>0 && y->count>0) || (xy->count>0))) {
      printf("\n%s\n\n", "User has to provide either x and y parameter or the xy_file parameter.");
      usage(argv[0],argtable);
    }

    if((xy->count>0) && !ezgdal_file_exists((char *)(xy->sval[0]))) {
      printf("\nFile [%s] does not exist.\n\n", xy->sval[0]);
      usage(argv[0],argtable);
    }

    dh_f = sml_open_layer((char *)(inp->sval[0]));
    if(dh_f==NULL) {
      printf("\nFile [%s] cannot be opened.\n\n", inp->sval[0]);
      usage(argv[0],argtable);
    }

    xybuf=(char*)calloc(1,dh_f->cell_size);
    if(app->count>0)
      f = fopen(out->sval[0],"a");
    else
      f = fopen(out->sval[0],"w");

    if(x->count>0 && y->count>0) {
      if(desc->count>0)
        description = (char *)(desc->sval[0]);

      sml_read_cell_from_layer_xy(dh_f,xybuf,x->dval[0],y->dval[0]);
      sml_write_cell_txt(f,x->dval[0],y->dval[0],description,xybuf,dh_f);
    } else {
      fxy = fopen((char *)xy->sval[0],"r");
      int line = 1;
      while(read_xy_txt(line++,fxy,&coord_x,&coord_y,desc_text, MAX_DESC_LEN)) {
        sml_read_cell_from_layer_xy(dh_f,xybuf,coord_x,coord_y);
        sml_write_cell_txt(f,coord_x,coord_y,desc_text,xybuf,dh_f);
      }
      fclose(fxy);
    }

    fclose(f);

    free(xybuf);
    sml_close_layer(dh_f);

printf("OK\n");

return 0;

}
