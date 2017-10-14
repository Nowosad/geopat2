#ifndef _TOOLS_COMMON_H_
#define _TOOLS_COMMON_H_

#include <stdio.h>

#define MAX_DESC_LEN 10240

/*
int file_exists(const char *fname);
void show_progress(int col, int cols);
void show_message(FILE *f,char *message);
*/
char *build_file_name(char *name, int idx);

#endif

