/****************************************************************************
 *
 * MODULE:	Signature normalization library
 * AUTHOR(S):	Pawel Netzel
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <stdlib.h>
#include <string.h>

#include "methods.h"
#include "methods_list.h"

normalization_func *get_normalization_method(char *method_name) {

  normalization_rec *p = normalizations_list;

  while(p->name != NULL) {
    if(strcmp(method_name,p->name)==0)
      return p->dist;
    p++;
  }
  
  return NULL;
}

char *get_normalization_description(char *method_name) {
  
  normalization_rec *p = normalizations_list;
  
  while(p->name != NULL) {
    if(strcmp(method_name,p->name)==0)
      return p->description;
    p++;
  }
  
  return NULL;
}

char *list_all_normalization_methods() {
  int len;
  char *buf;
  
  normalization_rec *p = normalizations_list;
  len = 0;

  while(p->name != NULL) {
    len += (int)strlen(p->name) + (int)strlen(p->description) + 5;
    p++;
  }
  buf = malloc(len);
  buf[0]='\0';
  p = normalizations_list;
  while(p->name != NULL) {
    strcat(buf,p->name);
    strcat(buf,"\t");
    strcat(buf,p->description);
    strcat(buf,"\n");
    p++;
  }

  return buf;

}