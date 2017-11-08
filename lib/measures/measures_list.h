#ifndef _MEASURES_LIST_H_
#define _MEASURES_LIST_H_

/****************************************************************************
*
* MODULE:	Similarity measures library (measures list)
* AUTHOR(S):	Pawel Netzel, Jakub Nowosad
* COPYRIGHT:	(C) Space Informatics Lab, Univeristy of Cincinnati
*
*		This program is free software under the GNU General Public
*		License (>=v2). Read the file COPYING that comes with GRASS
*		for details.
*
*****************************************************************************/

#include "measures.h"

/**
* distances of N vectors
*/
extern distance_func jensen_shannon;
extern distance_func wave_hedges;
extern distance_func triangular;

/**
* distances of 2 vectors
*/
extern distance_func euclidean;
extern distance_func euclidean_norm;
extern distance_func euclidean_period;
extern distance_func jaccard;
extern distance_func cosine;
extern distance_func rozicka;
extern distance_func rozickap;
extern distance_func hassanat;
extern distance_func ardiff;

/**
* distances of 2 time series
*/
extern distance_func tsEUC;
extern distance_func tsEUCP;
extern distance_func tsDTW;
extern distance_func tsDTWP;
extern distance_func tsDTWPa;

typedef struct {
        char *name;
        distance_func *dist;
        char *description;
} measure_rec;

measure_rec measures_list[] = {
        { "jsd",  jensen_shannon, "Jensen Shannon Divergence" },
        { "tri",  triangular, "Triangular"},
        { "euc",  euclidean,      "Euclidean distance" },
        { "eucn", euclidean_norm, "Normalized euclidean distance" },
        { "wh",   wave_hedges,    "Wave-Hedges distance" },
        { "jac",  jaccard,        "Jaccard distance" },
        /*************************
         *
         *   Experimental code
         * 
         { "eucp", euclidean_period, "Normalized euclidean distance (periodic)" },
         { "cos",  cosine,         "Cosine distance" },
         { "roz",  rozicka,        "Rozicka distance" },
         { "rozp", rozickap,       "Rozicka distance (extended)" },
         { "hass", hassanat,       "Hassanat distance" },
         *
         ************************/        
        { "tsEUC",   tsEUC,       "time series - euclidean distance" },
        { "tsEUCP",  tsEUCP,      "periodic time series - euclidean distance" },
        { "tsDTW",   tsDTW,       "time series - Dynamic Time Warping distance" },
        { "tsDTWP",  tsDTWP,      "time series - Periodic Dynamic Time Warping distance" },
        { "tsDTWPa", tsDTWPa,     "time series - Synchoronized Dynamic Time Warping distance" },
        { NULL, NULL, NULL }
};

#endif
