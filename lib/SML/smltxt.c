/****************************************************************************
 *
 * LIBRARY:	SML - Spatial Matrix Library
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	library for storing and accessing spatial grid of
 *		numbers, vectors, and tensors
 * COPYRIGHT:	(C) Pawel Netzel
 *              http://pawel.netzel.pl
 *
 *		This program is free software under the GNU Lesser General 
 *		Public License (>=v3). Read the file lgpl-3.0.txt
 *
 *****************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#define DLL_EXPORT

#include "sml.h"


void sml_write_cell_txt(FILE *f, double x, double y, char *desc, void *cell, SML_DATA_HEADER *dh) {

  int i;

  fprintf(f,"[%.10lf,%.10lf] \"%s\" %d:(%d",x,y,desc,dh->cell_type->dim,dh->cell_type->dims[0]);
  for(i=1; i<dh->cell_type->dim; i++)
    fprintf(f,",%d",dh->cell_type->dims[i]);
  if(!sml_is_cell_null(cell)) {
    fprintf(f,") => ");
    fprintf(f,"%.18lf",sml_get_cell_val_dbl(dh,cell,0));
    for(i=1; i<dh->cell_N_elements; i++)
      fprintf(f,",%.18lf",sml_get_cell_val_dbl(dh,cell,i));
    fprintf(f,"\n");
  } else {
    fprintf(f,") -> ");
    fprintf(f,"NO DATA\n");
  }
}

void sml_write_cell_csv(FILE *f, double x, double y, char *desc, void *cell, int use_nodata, double nodata, int decimals, SML_DATA_HEADER *dh) {

  int i;
  char format[10];

  if(sml_is_cell_null(cell) && !use_nodata) return;

  sprintf(format,",%%.%dlf",decimals);

  if(desc!=NULL)
    fprintf(f,"%s,",desc);
  fprintf(f,"%.10lf,%.10lf",x,y);
  if(!sml_is_cell_null(cell)) {
    for(i=0; i<dh->cell_N_elements; i++)
      fprintf(f,format,sml_get_cell_val_dbl(dh,cell,i));
  } else {
    for(i=0; i<dh->cell_N_elements; i++) 
      fprintf(f,format,nodata);
  }
  fprintf(f,"\n");
}

void sml_write_dblbuf_txt(FILE *f, double x, double y, char *desc, double *buffer, int size, int dim, int *dims) {

  int i;
  
  fprintf(f,"[%.10lf,%.10lf] \"%s\" %d:(%d",x,y,desc,dim,dims[0]);
  for(i=1; i<dim; i++)
    fprintf(f,",%d",dims[i]);
  if(size>0) {
    fprintf(f,") => ");
    fprintf(f,"%.18lf",buffer[0]);
    for(i=1; i<size; i++)
      fprintf(f,",%.18lf",buffer[i]);
    fprintf(f,"\n");
  } else {
    fprintf(f,") -> ");
    fprintf(f,"NO DATA\n");
  }
}

int sml_read_dblbuf_txt(FILE *f, double *x, double *y, char *desc, double *buffer, int size, SML_CELL_TYPE *ct) {

  int i = 0;
  double *d;
  char c;
  int s;

  desc[0] = '\0';
  i=fscanf(f,"[%lf,%lf] \"%[^\"]\" %d:(",x,y,desc,&(ct->dim));
  if(i!=4) {
    fscanf(f,"%*[^\n]");
    fscanf(f,"%*c");
    return -4;
  }

  if(ct->dims==NULL)
    ct->dims = (int *)malloc(ct->dim*sizeof(int));
  if(ct->len==NULL)
    ct->len = (int *)malloc(ct->dim*sizeof(int));
  ct->d_type = SML_DOUBLE;

  s = 1;
  for(i=0; i<ct->dim; i++) {
    fscanf(f,"%d%*c",&(ct->dims[i]));
    s *= ct->dims[i];
  }

/*
  ct->len[ct->dim-1] = 1;
  for(i=ct->dim-2; i>=0; i--)
    ct->len[i]=ct->len[i+1]*ct->dims[i];
*/


  if(size<=0)
    size = s;
  else if(size != s) {
    fscanf(f,"%*[^\n]");
    fscanf(f,"%*c");
    return -3;
  }


  s = size;
  fscanf(f,"%*c%c%*[^-1234567890.\n]",&c);

  if(c=='-') {
    fscanf(f,"%*c");
    return 0;
  }

  d = buffer;
  i = 0;
  while(i<size && fscanf(f,"%lf%c",d,&c)==2) {
    i++;
    if(c=='\n') {
      if(i<size) {
        fscanf(f,"%*[^\n]");
        fscanf(f,"%*c");
        return -1;
      }
    };
    d++;
  }
  if(c!='\n') {
    fscanf(f,"%*[^\n]");
    fscanf(f,"%*c");
    return -2;
  }

  return i;
}
