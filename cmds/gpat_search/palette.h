#ifndef _PALETTE_H_
#define _PALETTE_H_

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

typedef struct {
  int ncolors;
  double A, B;
  double (*pal)[5];
} PALETTE;

PALETTE *read_palette(const char *file_name);
PALETTE *create_default_palette();
void free_palette(PALETTE *p);

#endif