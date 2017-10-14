#ifndef _COMPATIBILITY_H_
#define _COMPATIBILITY_H_

/****************************************************************************
 *
 * MODULE:	Compatibility with GRASS GeoPAT
 * AUTHOR(S):	Pawel Netzel
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include <stdio.h>

#define NCINDEX(r,c,nc) (r)*(nc)+(c)
#define _(ARG) (ARG)

void G_message(const char *msg, ...);
void G_fatal_error(const char *msg, ...);
void G_warning(const char *msg, ...);

#endif