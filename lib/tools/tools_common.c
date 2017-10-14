#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
int file_exists(const char *fname) {
  FILE *f = fopen(fname,"r");
  if(f) {
    fclose(f);
    return 1;
  }
  return 0;
}

void show_progress(int col, int cols) {
//  if(fl_quiet) return;
  static int p = -1;
  int n = (int)(100.0*col/cols);
  if(p!=n) {
    p=n;
    printf("\b\b\b\b%3d%%",p); fflush(stdout);
    if(p==100) printf("\n");
  }
}

void show_message(FILE *f,char *message) {
//  if(fl_quiet) return;
  fprintf(f,"%s\n",message);
  fflush(f);
}
*/
char *build_file_name(char *name, int idx) {
  char *fname;
  if(strlen(name)==0) {
    fname = malloc(30*sizeof(char));
    sprintf(fname,"file_%d.tif",idx);
  } else {
    fname = malloc((strlen(name)+5)*sizeof(char));
    sprintf(fname,"%s.tif",name);
  }
  return fname;
}
