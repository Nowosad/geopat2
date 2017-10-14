#ifndef _SIGNATURES_LIST_H_
#define _SIGNATURES_LIST_H_

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
#include "signatures.h"

extern signature_func cartesianproduct;
extern signature_len_func cartesianproduct_len;
extern signature_func coocurrence;
extern signature_len_func coocurrence_len;
extern signature_func decomposition;
extern signature_len_func decomposition_len;
extern signature_func full_decomposition;
extern signature_len_func full_decomposition_len;
extern signature_func local_binary_pattern;
extern signature_len_func local_binary_pattern_len;
extern signature_func landind;
extern signature_len_func landind_len;
extern signature_func landind_short;
extern signature_len_func landind_short_len;
extern signature_func jcov;
extern signature_len_func jcov_len;
extern signature_func H;
extern signature_len_func H_len;

typedef struct {
	char *name;
	signature_func *signature;
	signature_len_func *signature_len;
	char *description;
} signature_rec;

signature_rec signatures_list[] = {
	{ "prod", cartesianproduct, cartesianproduct_len, "Cartesian product of input category lists" },
	{ "cooc", coocurrence, coocurrence_len, "Spatial coocurrence of categories" },
	{ "sdec", decomposition, decomposition_len, "Simple 2-level decomposition" },
	{ "fdec", full_decomposition, full_decomposition_len, "Full decomposition" },
	{ "lbp", local_binary_pattern, local_binary_pattern_len, "Histogram of local binary patterns" },
	{ "lind", landind, landind_len, "Landscape indices vector" },
	{ "linds", landind_short, landind_short_len, "Selected landscape indices vector" },
	{ "jcov", jcov, jcov_len, "J-Coocurrence vector" },
	{ "ent", H, H_len, "Shannon entropy" },
	{ NULL, NULL, NULL, 0 }
};

#endif