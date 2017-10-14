#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sml.h>

/*
void write_cell_txt(FILE *f, double x, double y, char *desc, void *cell, SML_DATA_HEADER *dh) {

  int i;
  
  fprintf(f,"[%.10lf,%.10lf] (%s) => ",x,y,desc);
  if(!sml_is_cell_null(cell)) {
    fprintf(f,"%.18lf",sml_get_cell_val_dbl(dh,cell,0));
    for(i=1; i<dh->cell_N_elements; i++)
      fprintf(f,",%.18lf",sml_get_cell_val_dbl(dh,cell,i));
    fprintf(f,"\n");
  } else
    fprintf(f,"NO DATA\n");

}

void write_dblbuf_txt(FILE *f, double x, double y, char *desc, double *buffer, int size) {

  int i;
  
  fprintf(f,"[%.10lf,%.10lf] (%s) => ",x,y,desc);
  if(size>0) {
    fprintf(f,",%.18lf",buffer[0]);
    for(i=1; i<size; i++)
      fprintf(f,",%.18lf",buffer[i]);
    fprintf(f,"\n");
  } else
    fprintf(f,"NO DATA\n");

}

int read_dblbuf_txt(FILE *f, double *x, double *y, char *desc, double *buffer, int size) {

  int i = 0;
  double *d;
  char c;
  
  fscanf(f,"[%lf,%lf] (%[^)]",x,y,desc);

//printf("x: %lf, y: %lf, desc: %s\n",*x,*y,desc); fflush(stdout);

  fscanf(f,"%*[^-1234567890.\n]");
  d = buffer;
  while(i<size && fscanf(f,"%lf%c",d,&c)==2) {
//printf("v(%d): %lf\n",i,*d); fflush(stdout);
    i++;
    if(c=='\n') break;
    d++;
  }

  return i;

}
*/
/*
int read_xy_txt(FILE *f, double *x, double *y, char *desc, int max_size) {
  char format[20];
  int i = fscanf(f,"%lf,%lf",x,y);
  int n;
  max_size--;
  sprintf(format,"%%%d[^,]%%n",max_size);
  if(i!=2) {
    i = fscanf(f,format,desc,&n);
    if(n==max_size) fscanf(f,"%*[^,]");
    fscanf(f,"%*c");
    if(i == 1)
      i = fscanf(f,"%lf,%lf",x,y);
  } else
    *desc = '\0';
  fscanf(f,"%*[^\n]");
  fscanf(f,"%*c");
  return i == 2;

}
*/

char *create_fname(char *desc) {
  if(desc==NULL) return NULL;
  int i, n;
  char c, *res;
  n = (int)strlen(desc);
  if(n==0) return NULL;
  for(i=0; i<n; i++) {
    c = desc[i];
    if(!isalnum(c) && !(c=='.' || c=='-' || c=='_' || c=='/' || c=='\\' || c==':'))
      desc[i]='_';
  }
  res = (char *)malloc((n+6)*sizeof(char));
  sprintf(res,"%s.tif",desc);
  return res;
}

int read_xy_txt(int line, FILE *f, double *x, double *y, char *desc, int max_size) {
  char format[20];
  int i = fscanf(f,"%lf,%lf",x,y);
  desc[0] = '\0';
  sprintf(format,",%%%d[^\n]%%n",max_size-1);
  fscanf(f,format,desc);
  if(strlen(desc)==0)
    sprintf(desc,"loc %05d",line);
  fscanf(f,"%*c");
  return i == 2;
}
