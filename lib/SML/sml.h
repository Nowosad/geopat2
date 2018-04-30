#ifndef _SML_H_
#define _SML_H_

/****************************************************************************
 *
 * LIBRARY:	SML - Spatial Matrix Library
 * AUTHOR(S):	Pawel Netzel
 * PURPOSE:	library for storing and accessing spatial grid of
 *		numbers, vectors, and tensors
 * COPYRIGHT:	(C) Pawel Netzel
 *              http://pawel.netzel.pl
 *
 *		This program is free software under the GNU Lesser General 
 *		Public License (>=v3). Read the file lgpl-3.0.txt
 *
 *****************************************************************************/



#include <stdio.h>
#include <stdlib.h>

#ifdef _MSC_VER
  #ifdef DLL_EXPORT
    #define SML_DLL_API __declspec(dllexport)
  #else 
    #define SML_DLL_API __declspec(dllimport)
  #endif
#else
  #define SML_DLL_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define SML_MAX_BUF 128
#define SML_DATA_NAME_LEN 4096
#define SML_DATA_DESC_LEN 32768
#define SML_DATA_PROJ_LEN 8192

typedef enum SML_FILE_STATUS {
		SML_NEW,
		SML_EXISTING
	} SML_FILE_STATUS;


typedef enum SML_D_TYPE {
		SML_BYTE,
		SML_SHORT_INT,
		SML_INT,
		SML_FLOAT,
		SML_DOUBLE
	} SML_D_TYPE;

typedef double SML_AFFINE_TRANSFORM[6];

typedef struct SML_WINDOW {
	SML_AFFINE_TRANSFORM at;
	int rows;
	int cols;
	char proj[SML_DATA_PROJ_LEN];
} SML_WINDOW;

typedef struct SML_CELL_TYPE {
	int dim;
	int *dims;
	int *len;
	SML_D_TYPE d_type;
} SML_CELL_TYPE;

typedef struct SML_DATA_HEADER {
	FILE *f;
	SML_FILE_STATUS file_status;
	char name[SML_DATA_NAME_LEN];
	char desc[SML_DATA_DESC_LEN];
	SML_CELL_TYPE *cell_type;
	SML_WINDOW *file_win;
	int cell_N_elements;
	int cell_element_size;
	int cell_size;
} SML_DATA_HEADER;


/**   WINDOW
 * 
 */
SML_DLL_API SML_WINDOW *sml_create_window(int rows, int cols, double at0, double at1, double at2, double at3, double at4, double at5, char* proj);
SML_DLL_API SML_WINDOW *sml_create_window_copy(SML_WINDOW *w);
SML_DLL_API void sml_free_window(SML_WINDOW *w);

/**   CELL_TYPE
 * 
 */
SML_DLL_API SML_CELL_TYPE *sml_create_cell_type(SML_D_TYPE d, int ndim, int* dims);
SML_DLL_API SML_CELL_TYPE *sml_create_cell_type_copy(SML_CELL_TYPE *ct0);
SML_DLL_API void sml_free_cell_type(SML_CELL_TYPE *ct);

/**   LAYER
 * 
 */
SML_DLL_API SML_DATA_HEADER* sml_open_layer(char* fname);
SML_DLL_API SML_DATA_HEADER* sml_create_layer(
		char* fname,
		SML_CELL_TYPE *cell_type,
		SML_WINDOW *w);
SML_DLL_API void sml_close_layer(SML_DATA_HEADER *dh);
SML_DLL_API void sml_set_layer_description(SML_DATA_HEADER *dh, char *desc[], int cnt);
SML_DLL_API void sml_read_cell_from_layer(SML_DATA_HEADER *dh, void *cell, int col, int row);
SML_DLL_API void sml_read_cell_from_layer_xy(SML_DATA_HEADER *dh, void *cell, double x, double y);
SML_DLL_API void sml_read_row_from_layer(SML_DATA_HEADER *dh, void *cell_row, int row);
SML_DLL_API void sml_write_row_to_layer(SML_DATA_HEADER *dh, void *cell_row, int row);
SML_DLL_API void sml_write_next_row_to_layer(SML_DATA_HEADER *dh, void *cell_row);

/**   CELL BUFFER (row)
 * 
 */
SML_DLL_API void *sml_create_cell_row_buffer(SML_DATA_HEADER *dh);
SML_DLL_API void *sml_get_cell_pointer(SML_DATA_HEADER *dh, void *row_buffer, int i);

/**   CELL
 * 
 */
SML_DLL_API char sml_is_cell_null(void *cell);
SML_DLL_API void sml_set_cell_null(void *cell);
SML_DLL_API void sml_set_cell_not_null(void *cell);
SML_DLL_API void *sml_get_cell_data(void *cell);
SML_DLL_API double sml_get_cell_val_dbl(SML_DATA_HEADER *dh, void *cell, int i);
SML_DLL_API int sml_get_cell_val_int(SML_DATA_HEADER *dh, void *cell, int i);
SML_DLL_API void *sml_get_cell_val(SML_DATA_HEADER *dh, void *cell, int i);
SML_DLL_API void sml_set_cell_val_dbl(SML_DATA_HEADER *dh, double v, void *cell, int i);
SML_DLL_API void sml_set_cell_val_int(SML_DATA_HEADER *dh, int v, void *cell, int i);

/**   CELL - text file
 * 
 */
SML_DLL_API void sml_write_cell_txt(FILE *f, double x, double y, char *desc, void *cell, SML_DATA_HEADER *dh);
SML_DLL_API void sml_write_dblbuf_txt(FILE *f, double x, double y, char *desc, double *buffer, int size, int dim, int *dims);
SML_DLL_API void sml_write_cell_csv(FILE *f, double x, double y, char *desc, void *cell, int use_nodata, double nodata, int decimals, SML_DATA_HEADER *dh);
SML_DLL_API int sml_read_dblbuf_txt(FILE *f, double *x, double *y, char *desc, double *buffer, int size, SML_CELL_TYPE *ct);


/**   COORDS -> ROW,COL
 * 
 */
SML_DLL_API int sml_xy2c(SML_DATA_HEADER *dh, double x, double y);
SML_DLL_API int sml_xy2r(SML_DATA_HEADER *dh, double x, double y);


/**   COL,ROW -> COORDS
 * 
 */
SML_DLL_API double sml_cr2y(SML_DATA_HEADER *dh, int col, int row);
SML_DLL_API double sml_cr2x(SML_DATA_HEADER *dh, int col, int row);

#ifdef __cplusplus
}
#endif


#endif
