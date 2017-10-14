/****************************************************************************
 *
 * PROGRAM:	gpat_search - part of GeoPAT 2
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	program for calculating similarity layer;
 *		functionality based on p.sim.search from
 *		GRASS GeoPAT by Jasiewicz, Netzel, Stepinski
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
#include <assert.h>

#include "palette.h"

/**************************************
 *
 *  reading color palette from file
 *
 *  file format:
 *
 *    PALETTE
 *    ncolors scaleA scaleB
 *    comment (if no commnts, leave empty line)
 *    idx0 r0 g0 b0
 *    idx1 r1 g1 b1
 *    idx2 r2 g2 b2
 *    ........
 *
 ***************************************/

PALETTE *read_palette(const char *file_name) {
  PALETTE *p = NULL;
  int i,n,k;
  FILE *f;
  char buf[1024];
  double a,b;

  f = fopen(file_name,"r");
  assert(f!=NULL);

  fscanf(f,"%10[^\n]%*c",buf);
  if(strcmp(buf,"PALETTE")==0) {

    fscanf(f,"%d %lf %lf\n", &n, &a, &b);
    assert(n>0 && n<256);
    assert(a!=0);

    fscanf(f,"%*[^\n]%*c");

    p = malloc(sizeof(PALETTE));
    assert(p!=NULL);

    p->ncolors = n;
    p->A = a;
    p->B = b;
    p->pal = (double (*)[5])malloc(n*sizeof(p->pal[0]));
    assert(p->pal!=NULL);
/*
    for(i=0; i<n; i++)
      p->pal[i] = malloc(4*sizeof(double));
*/
    for(i=0; i<n; i++) {
      k = fscanf(f,"%lf %lf %lf %lf %lf",&(p->pal[i][0]),&(p->pal[i][1]),&(p->pal[i][2]),&(p->pal[i][3]),&(p->pal[i][4]));
      assert(k==5);
    }

  }

  fclose(f);

  return p;
}

PALETTE *create_default_palette() {
  PALETTE *p = NULL;
  int i;

  p = malloc(sizeof(PALETTE));
  assert(p!=NULL);

  p->ncolors = 256;
  p->A = 255.0;
  p->B = 0.0;
  p->pal = (double (*)[5])malloc(256*sizeof(p->pal[0]));
  assert(p->pal!=NULL);
/*
  for(i=0; i<n; i++) {
    p->pal[i] = malloc(4*sizeof(double));
    assert(p->pal[i]!=NULL);
  }
*/
  for(i=0; i<256; i++) {
    p->pal[i][0] = i;     // element no / value
    p->pal[i][1] = 255-i; // R - (0, 255)
    p->pal[i][2] = i;     // G - (0, 255)
    p->pal[i][3] = 0;     // B - (0, 255)
    p->pal[i][4] = 255;   // opacity - (0, 255)
  }

  return p;
}

void free_palette(PALETTE *p) {
  if(p==NULL) return;
//  int i;
//  for(i=0; i<p->ncolors; i++)
  free(p->pal);
  free(p);
}

