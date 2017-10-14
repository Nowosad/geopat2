#ifndef _MEASURES_H_
#define _MEASURES_H_

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
#include <stdarg.h>

typedef double distance_func(double**, int, int, int, int*, ...);

distance_func *get_distance(char *distance_name);
char *get_distance_description(char *distance_name);
char *list_all_distances();

#endif
