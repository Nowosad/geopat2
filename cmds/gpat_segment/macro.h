#ifndef __MACRO_H__
#define __MACRO_H__

/****************************************************************************
 *
 * MODULE:	segmentation of motifels grid
 * AUTHOR(S):	Jaroslaw Jasiewicz, Jacek Niesterowicz, Tomasz Stepinski
 * PURPOSE:	information retrival using categorical maps:
 *		compares grid of histograms
 * COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
 *
 *		This program is free software under the GNU General Public
 *		License (>=v2). Read the file COPYING that comes with GRASS
 *		for details.
 *
 *****************************************************************************/

#define NR(x) r + nextr[(x)]
#define NRR(x) r + row + nextr[(x)]
#define NC(x) c + nextc[(x)]
#define NCC(x) c+col + nextc[(x)]
/*
static int nextr[] = {0, 1, 1,  1};
static int nextc[] = {1, 1, 0, -1};

static int Snextr[9] = { 0, -1, -1, -1, 0, 1, 1, 1, 0 };
static int Snextc[9] = { 0, 1, 0, -1, -1, -1, 0, 1, 1 };
*/
#define SNR(x) r + Snextr[(x)]
#define SNC(x) c + Snextc[(x)]
//#define SNCC(x) c+col + Snextc[(x)]
#define INDEX(r,c) (r)*ncols+(c)
#define NCINDEX(r,c,nc) (r)*(nc)+(c)
#define NOT_IN_REGION(x) (r+Snextr[(x)] < 0 || r+Snextr[(x)] > (nrows-1) || c+Snextc[(x)] < 0 || c+Snextc[(x)] > (ncols-1))

#define INF 1E+30
#endif
