/****************************************************************************
 *
 * PROGRAM:	gpat_grd2txt - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	program for converting binary grid of motifels to text format;
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
#include <assert.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "../../lib/ezGDAL/ezgdal.h"
#include "../../lib/SML/sml.h"

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
    SML_DATA_HEADER *dh;
    void *rowbuf;

    double x,y;
    int r,c;
    char desc[1024];
    void *cell;

    FILE *f;

    struct arg_str  *inp   = arg_str1("i","input","<file_name>","name of input file (GRID)");
    struct arg_str  *out   = arg_str1("o","output","<file_name>","name of output file (TXT)");
    struct arg_lit  *help  = arg_lit0("h","help","print this help and exit");
    struct arg_end  *end   = arg_end(20);
    void* argtable[] = {inp,out,help,end};

    int nerrors = arg_parse(argc,argv,argtable);

    if (help->count > 0) 
      usage(argv[0],argtable);

    /* If the parser returned any errors then display them and exit */
    if (nerrors > 0) {
      /* Display the error details contained in the arg_end struct.*/
      arg_print_errors(stdout,end,argv[0]);
      usage(argv[0],argtable);
    }

    if(inp->count>0 && !ezgdal_file_exists((char *)(inp->sval[0]))) {
      printf("\nFile '%s' does not exists!\n\n", inp->sval[0]);
      usage(argv[0],argtable);
    }

    dh = sml_open_layer((char *)(inp->sval[0]));
    rowbuf = sml_create_cell_row_buffer(dh);

    f = fopen(out->sval[0],"w");

    ezgdal_show_progress(stdout,0,dh->file_win->rows);
    for(r=0; r<dh->file_win->rows; r++) {

      ezgdal_show_progress(stdout,r,dh->file_win->rows);
      sml_read_row_from_layer(dh,rowbuf,r);

      for(c=0; c<dh->file_win->cols; c++) {
        cell = sml_get_cell_pointer(dh,rowbuf,c);
        if(!sml_is_cell_null(cell)) {
          x = sml_cr2x(dh,c,r);
          y = sml_cr2y(dh,c,r);
          sprintf(desc,"%s_%d_%d",(char *)(inp->sval[0]),c,r);
          sml_write_cell_txt(f, x, y, desc, cell, dh);
        }
      }
      fflush(f);
    }
    ezgdal_show_progress(stdout,100,100);

    fclose(f);

    free(rowbuf);
    sml_close_layer(dh);

return 0;

}

