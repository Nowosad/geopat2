#ifndef _TOOLS_CALL_TXT_H_
#define _TOOLS_CALL_TXT_H_

#include <stdio.h>
#include "../../lib/SML/sml.h"

/*
void write_cell_txt(FILE *f, double x, double y, char *desc, void *cell, SML_DATA_HEADER *dh);
void write_dblbuf_txt(FILE *f, double x, double y, char *desc, double *buffer, int size);
int read_dblbuf_txt(FILE *f, double *x, double *y, char *desc, double *buffer, int size);
*/
int read_xy_txt(int line, FILE *f, double *x, double *y, char *desc, int max_size);
char *create_fname(char *desc);

#endif
