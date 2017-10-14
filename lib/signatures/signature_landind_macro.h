#ifndef __MACRO_H__
#define __MACRO_H__

/****************************************************************************
 *
 * MODULE:	landscape indices signature
 * AUTHOR(S):	Jacek Niesterowicz
 * PURPOSE:	calculating a vector of landscape indices:
 *
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/

#define NR(x) r + nextr[(x)]
#define NC(x) c + nextc[(x)]
#define NCC(x) c+col + nextc[(x)]
/*
static long int nextr[] = {0, 1, 1,  1};
static long int nextc[] = {1, 1, 0, -1};
*/
static long int Snextr[9] = { 0, -1, -1, -1, 0, 1, 1, 1, 0 };
static long int Snextc[9] = { 0, 1, 0, -1, -1, -1, 0, 1, 1 };

static int four_conn[9] = {0, 0, 1, 0, 1, 0, 1, 0, 1};

#define SNR(x) r + Snextr[(x)]
#define SNC(x) c + Snextc[(x)]
#define INDEX(r,c) (r)*ncols+(c)
#define NCINDEX(r,c,nc) (r)*(nc)+(c)
#define NOT_IN_REGION(x) (r+Snextr[(x)] < 0 || r+Snextr[(x)] > (nrows-1) || c+Snextc[(x)] < 0 || c+Snextc[(x)] > (ncols-1))
#define IS_FOUR_CONN(x) four_conn[(x)]
#define INF 1E+30

#endif
