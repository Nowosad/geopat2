#ifndef _SIGNATURES_H_
#define _SIGNATURES_H_

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

#include <stdarg.h>
#include <ezgdal.h>

typedef int signature_func(EZGDAL_FRAME**, int, double*, int, ...);
typedef int signature_len_func(EZGDAL_LAYER**, int, ...);

signature_func *get_signature(char *signature_name);
signature_len_func *get_signature_len(char *signature_name);
char *get_signature_description(char *signature_name);
char *list_all_signatures();

#endif
