#include <stdio.h>
#include "palette.h"

void main() {

  int i;

  PALETTE *p = read_palette("palette.txt");
//  PALETTE *p = create_default_palette();

  printf("%d %lf %lf\n",p->ncolors, p->A, p->B);
  for(i=0; i<p->ncolors; i++) 
    printf("  %d - %lf %lf %lf %lf\n",i,p->pal[i][0],p->pal[i][1],p->pal[i][2],p->pal[i][3]);

  free_palette(p);
}