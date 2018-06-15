#ifndef _EZGDAL_H_
#define _EZGDAL_H_

/****************************************************************************
 *
 * LIBRARY:	ezGDAL - easy GDAL access
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	library for reading GDAL layers, writting GeoTIFFs, processing
 *		data row-by-row, stripe-by-stripe, frame-by-frame
 * COPYRIGHT:	(C) Pawel Netzel
 *              http://pawel.netzel.pl
 *
 *		This program is free software under the GNU Lesser General 
 *		Public License (>=v3). Read the file lgpl-3.0.txt
 *
 *****************************************************************************/

#include <stdio.h>
#include <gdal.h>

#ifdef _MSC_VER
  #ifdef DLL_EXPORT
    #define EZGDAL_DLL_API __declspec(dllexport)
  #else 
    #define EZGDAL_DLL_API __declspec(dllimport)
  #endif
#else
  #define EZGDAL_DLL_API
#endif

#ifdef __cplusplus
extern "C" {
#endif


#define EZGDAL_FRAMESET_LEN 10240
#define EZGDAL_FRAMESET_STEP 10240
#define MAX_FRAME_BUFFER_SIZE 4294967295

typedef struct {
  double min,max,avg,std;
  double hist_min,hist_max,hist_step;
  int hist_N;
  int map_max_val;
  int *map_cat;
  GUIntBig *hist;
} EZGDAL_STATS;

typedef struct EZGDAL_STRIPE EZGDAL_STRIPE;
typedef struct EZGDAL_FRAMESET EZGDAL_FRAMESET;
typedef struct EZGDAL_LAYER EZGDAL_LAYER;

union EZGDAL_FRAME_OWNER {
  EZGDAL_STRIPE *stripe;
  EZGDAL_FRAMESET *frameset;
};

typedef struct {
  union EZGDAL_FRAME_OWNER owner;
  int cols, rows;
  int col1, col2;
  int row1, row2;
  double *private_buffer;
  double **buffer;
} EZGDAL_FRAME;

struct EZGDAL_FRAMESET {
  EZGDAL_LAYER *layer;
  int frames;
  int frameset_len;
  EZGDAL_FRAME **frame;
};

struct EZGDAL_STRIPE {
  EZGDAL_LAYER *layer;
  int row1, row2;
  int rows;
  double **buffer;
  int frames;
  EZGDAL_FRAME *frame;
};

struct EZGDAL_LAYER {
  GDALDatasetH dataset_h;
  GDALRasterBandH band_h;
  int rows;
  int cols;
  int is_no_data;
  double no_data;
  double *buffer;
  EZGDAL_STRIPE *stripe;
  EZGDAL_FRAMESET *frameset;
  EZGDAL_STATS *stats;
};

/*==========================================*/
/*              TOOLS                       */

EZGDAL_DLL_API int  ezgdal_file_exists(const char *fname);
EZGDAL_DLL_API void  ezgdal_show_progress(FILE *f, int i, int N);
EZGDAL_DLL_API void  ezgdal_show_message(FILE *f, char *message);

EZGDAL_DLL_API int  ezgdal_is_bbox_ok(EZGDAL_LAYER **inputs, int ninputs);
EZGDAL_DLL_API int  ezgdal_is_projection_ok(EZGDAL_LAYER **inputs, int ninputs);

EZGDAL_DLL_API GDALDataType  ezgdal_data_type(char *data);

EZGDAL_DLL_API int  ezgdal_xy2c(EZGDAL_LAYER *layer, double x, double y);
EZGDAL_DLL_API double  ezgdal_cr2x(EZGDAL_LAYER *layer, int col, int row);
EZGDAL_DLL_API int  ezgdal_xy2r(EZGDAL_LAYER *layer, double x, double y);
EZGDAL_DLL_API double  ezgdal_cr2y(EZGDAL_LAYER *layer, int col, int row);

/*==========================================*/
/*              EZGDAL_LAYER                       */

EZGDAL_DLL_API EZGDAL_LAYER*  ezgdal_create_layer(char *fname, 
                         char *wkt,
                         char *d_type,
                         double *geo_tr,
                         int rows,
                         int cols,
                         int *nodata);
EZGDAL_DLL_API EZGDAL_LAYER*  ezgdal_open_layer(char *fname);
EZGDAL_DLL_API void  ezgdal_close_layer(EZGDAL_LAYER *layer);

EZGDAL_DLL_API void  ezgdal_set_palette255(EZGDAL_LAYER *layer, double palette[][5], int n);

EZGDAL_DLL_API void  ezgdal_read_buffer(EZGDAL_LAYER *layer, int row);
EZGDAL_DLL_API void  ezgdal_write_buffer(EZGDAL_LAYER *layer, int row);

EZGDAL_DLL_API char*  ezgdal_layer_get_wkt(EZGDAL_LAYER *layer);
EZGDAL_DLL_API double*  ezgdal_layer_get_at(EZGDAL_LAYER *layer);

EZGDAL_DLL_API int  ezgdal_is_null(EZGDAL_LAYER *layer, double v);
EZGDAL_DLL_API void  ezgdal_set_null(EZGDAL_LAYER *layer, double *v);

/*==========================================*/
/*              STATS                       */

EZGDAL_DLL_API void  ezgdal_calc_layer_stats(EZGDAL_LAYER *layer);
EZGDAL_DLL_API void  ezgdal_calc_value_map(EZGDAL_LAYER *layer, double min, double max, int N);
EZGDAL_DLL_API int  ezgdal_get_value_index(EZGDAL_LAYER *layer, double val);
EZGDAL_DLL_API double  ezgdal_get_index_value(EZGDAL_LAYER *layer, int idx);

/*==========================================*/
/*         STRIPE & FRAME                   */

EZGDAL_DLL_API EZGDAL_STRIPE*  ezgdal_create_stripe(EZGDAL_LAYER *layer, int row1, int row2);
EZGDAL_DLL_API EZGDAL_FRAME*  ezgdal_create_frame(EZGDAL_STRIPE *stripe, int col1);
EZGDAL_DLL_API int  ezgdal_create_all_frames(EZGDAL_STRIPE *stripe, int start, int shift);
EZGDAL_DLL_API void  ezgdal_free_stripe(EZGDAL_LAYER *layer);
EZGDAL_DLL_API void  ezgdal_free_all_frames(EZGDAL_STRIPE *stripe);

EZGDAL_DLL_API EZGDAL_FRAME*  ezgdal_get_frame(EZGDAL_STRIPE *stripe, int idx);
EZGDAL_DLL_API void  ezgdal_shift_frame_pos(EZGDAL_FRAME *frame, int col);

EZGDAL_DLL_API int  ezgdal_load_stripe_data(EZGDAL_STRIPE *stripe, int row1);
EZGDAL_DLL_API int  ezgdal_save_stripe_data(EZGDAL_STRIPE *stripe);

/*==========================================*/
/*         FRAMESET & FRAME                 */

EZGDAL_DLL_API void ezgdal_frameset_max_buffer_size(unsigned long size);
EZGDAL_DLL_API EZGDAL_FRAMESET*  ezgdal_create_frameset(EZGDAL_LAYER *layer);
EZGDAL_DLL_API EZGDAL_FRAMESET*  ezgdal_create_frameset_with_size(EZGDAL_LAYER *layer, int size);
EZGDAL_DLL_API EZGDAL_FRAME*  ezgdal_add_frameset_frame(EZGDAL_FRAMESET *frameset, int col1, int col2, int row1, int row2);
EZGDAL_DLL_API void  ezgdal_frameset_set_frame(EZGDAL_FRAME *frame, int col1, int col2, int row1, int row2);
EZGDAL_DLL_API void  ezgdal_free_frameset(EZGDAL_FRAMESET *frameset);
EZGDAL_DLL_API void  ezgdal_free_frameset_frame(EZGDAL_FRAME *frame);
EZGDAL_DLL_API void  ezgdal_free_frameset_all_frames(EZGDAL_FRAMESET *frameset);

EZGDAL_DLL_API EZGDAL_FRAME*  ezgdal_get_frameset_frame(EZGDAL_FRAMESET *frameset, int idx);
EZGDAL_DLL_API void  ezgdal_load_frameset_frame_data(EZGDAL_FRAME *frame);
EZGDAL_DLL_API void  ezgdal_unload_frameset_frame_data(EZGDAL_FRAME *frame);

#ifdef __cplusplus
}
#endif


#endif

