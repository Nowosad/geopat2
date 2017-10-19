#ifndef _METHODS_LIST_H_
#define _METHODS_LIST_H_

/****************************************************************************
 *
 * MODULE:	Signature normalization library
 * AUTHOR(S):	Pawel Netzel
 * COPYRIGHT:	(C) Space Informatics Lab, University of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/
#include "methods.h"

extern normalization_func normalization_01;
extern normalization_func normalization_pdf;
extern normalization_func normalization_N01;
extern normalization_func normalization_none;

typedef struct {
	char *name;
	normalization_func *dist;
	char *description;
} normalization_rec;

normalization_rec normalizations_list[] = {
	{ "01",   normalization_01,   "Normalize to the interval [0, 1]" },
	{ "pdf",  normalization_pdf,  "Normalize to pdf (sum(hi) = 1)" },
	{ "N01",  normalization_N01,  "Normalize to N(0, 1) (hi=(hi-avg)/std)" },
	{ "none", normalization_none, "Signature without any normalization" },
	{ NULL, NULL, NULL }
};

#endif
