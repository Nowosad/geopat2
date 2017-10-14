/****************************************************************************
 *
 * MODULE:	Similarity measures library
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

#include "measures.h"
#include "measures_list.h"

distance_func *get_distance(char *distance_name) {
  
  measure_rec *p = measures_list;
  
  while(p->name != NULL) {
    if(strcmp(distance_name,p->name)==0)
      return p->dist;
    p++;
  }
  
  return NULL;
}

char *get_distance_description(char *distance_name) {
  
  measure_rec *p = measures_list;
  
  while(p->name != NULL) {
    if(strcmp(distance_name,p->name)==0)
      return p->description;
    p++;
  }
  
  return NULL;
}

char *list_all_distances() {
  int len;
  char *buf;
  
  measure_rec *p = measures_list;
  len = 0;

  while(p->name != NULL) {
    len += (int)strlen(p->name) + (int)strlen(p->description) + 5;
    p++;
  }
  buf = malloc(len);
  buf[0]='\0';
  p = measures_list;
  while(p->name != NULL) {
    strcat(buf,p->name);
    strcat(buf,"\t");
    strcat(buf,p->description);
    strcat(buf,"\n");
    p++;
  }

  return buf;

}