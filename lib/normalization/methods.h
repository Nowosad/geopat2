#ifndef _METHODS_H_
#define _METHODS_H_

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

typedef int normalization_func(double*, int);

normalization_func *get_normalization_method(char *method_name);
char *get_normalization_description(char *method_name);
char *list_all_normalization_methods();

#endif
