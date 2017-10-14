/****************************************************************************
 *
 * MODULE:	Signatures library
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

#include "signatures.h"
#include "signatures_list.h"

signature_func *get_signature(char *signature_name) {
  
  signature_rec *p = signatures_list;
  
  while(p->name != NULL) {
    if(strcmp(signature_name,p->name)==0)
      return p->signature;
    p++;
  }
  
  return NULL;
}

signature_len_func *get_signature_len(char *signature_name) {
  
  signature_rec *p = signatures_list;
  
  while(p->name != NULL) {
    if(strcmp(signature_name,p->name)==0)
      return p->signature_len;
    p++;
  }
  
  return NULL;
}

char *get_signature_description(char *signature_name) {
  
  signature_rec *p = signatures_list;
  
  while(p->name != NULL) {
    if(strcmp(signature_name,p->name)==0)
      return p->description;
    p++;
  }
  
  return NULL;
}

char *list_all_signatures() {
  int len;
  char *buf;
  
  signature_rec *p = signatures_list;
  len = 0;

  while(p->name != NULL) {
    len += (int)strlen(p->name) + (int)strlen(p->description) + 5;
    p++;
  }
  buf = malloc(len);
  buf[0]='\0';
  p = signatures_list;
  while(p->name != NULL) {
    strcat(buf,p->name);
    strcat(buf,"\t");
    strcat(buf,p->description);
    strcat(buf,"\n");
    p++;
  }

  return buf;

}