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



#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <gdal.h>
#include <ogr_srs_api.h>
#include <cpl_string.h>


#define DLL_EXPORT

#include "ezgdal.h"


/*==========================================*/
/*              UTILS                       */
/*                                          */

static unsigned long max_frame_buffer_size = MAX_FRAME_BUFFER_SIZE;

int  ezgdal_file_exists(const char *fname) {
  FILE *f = fopen(fname,"r");
  if(f) {
    fclose(f);
    return 1;
  }
  return 0;
}

void  ezgdal_show_progress(FILE *f, int i, int N) {
  static int percent = -1;
  int n;
  if(N==0 || i<0 || i>N) 
	  n = 0;
  else 
	  n = (int)(100.0*(double)i/(double)N);
  if(percent!=n) {
    percent=n;
    fprintf(f, "\b\b\b\b%3d%%",percent);
	fflush(f);
    if(percent==100) 
		fprintf(f, "\n");
  }
}

void  ezgdal_show_message(FILE *f, char *message) {
  fprintf(f,"%s\n",message);
  fflush(f);
}


/*==========================================*/
/*                 TOOLS                    */
/*                                          */

int  ezgdal_is_bbox_ok(EZGDAL_LAYER **inputs, int ninputs) {
  double eps = 0.00000000001;
  double geo_transform[2][6];
  int i,j;
  
  if(ninputs<2) return 1;
  GDALGetGeoTransform(inputs[0]->dataset_h,geo_transform[0]);
  for(i=1; i<ninputs; i++) {
	GDALGetGeoTransform(inputs[i]->dataset_h,geo_transform[1]);
	for(j=0; j<6; j++)
	  if(fabs(geo_transform[0][j]-geo_transform[1][j])>eps)
	    return 0;
  }
  return 1;
}

int  ezgdal_is_projection_ok(EZGDAL_LAYER **inputs, int ninputs) {
  int i = 1;
  int proj_ok = TRUE;
  char *p, *p1;
  OGRSpatialReferenceH ref, ref1;

  if(ninputs<2) 
    return TRUE;

  p = (char *)GDALGetProjectionRef(inputs[0]->dataset_h);

  ref = OSRNewSpatialReference(p);

  while(proj_ok && i<ninputs) {
    p1 = (char *)GDALGetProjectionRef(inputs[i]->dataset_h);
    ref1 = OSRNewSpatialReference(p1);
    proj_ok = proj_ok && ((p==NULL && p1==NULL) || (p!=NULL && p1!=NULL && OSRIsSame(ref,ref1)));

    OSRDestroySpatialReference(ref1);
    i++;
  }

  OSRDestroySpatialReference(ref);
  return proj_ok;
}

GDALDataType  ezgdal_data_type(char *data) {
  if(data==NULL)
    return GDT_Float64;
  if(strcmp(data,"Byte")==0)
    return GDT_Byte;
  else if(strcmp(data,"UInt16")==0)
    return GDT_UInt16;
  else if(strcmp(data,"Int16")==0)
    return GDT_Int16;
  else if(strcmp(data,"UInt32")==0)
    return GDT_UInt32;
  else if(strcmp(data,"Int32")==0)
    return GDT_Int32;
  else if(strcmp(data,"Float32")==0)
    return GDT_Float32;
  else if(strcmp(data,"Float64")==0)
    return GDT_Float64;
  else
    return GDT_Float64;
}

char*  ezgdal_layer_get_wkt(EZGDAL_LAYER *layer) {
  char *p = (char *)GDALGetProjectionRef(layer->dataset_h);
  return p;  
}

double*  ezgdal_layer_get_at(EZGDAL_LAYER *layer) {
  double *p = malloc(6*sizeof(double));

  GDALGetGeoTransform(layer->dataset_h, p);
  return p;
}

int  ezgdal_xy2c(EZGDAL_LAYER *layer, double x, double y) {
	double *a = ezgdal_layer_get_at(layer);
        return (int)floor((x - a[0])/a[1]);


/*
	double c0 = a[2]*a[4]-a[1]*a[5];
	if(a[2]!=0.0 && c0!=0.0) 
	    return (int)((a[0]*a[5]-a[2]*a[3]+a[2]*y-a[5]*x)/c0);
	else if(a[2]==0 && a[1]!=0.0 && a[5]!=0.0)
	    return (int)((x-a[0])/a[1]);
	else 
	    return -1;
*/
};
double  ezgdal_cr2x(EZGDAL_LAYER *layer, int col, int row) { 
	double *a = ezgdal_layer_get_at(layer);
	return a[0]+a[1]*col+a[2]*row;
};

int  ezgdal_xy2r(EZGDAL_LAYER *layer, double x, double y) {
	double *a = ezgdal_layer_get_at(layer);
        return (int)floor((y - a[3])/a[5]);

/*
	double c0 = a[2]*a[4]-a[1]*a[5];
	if(a[2]!=0.0 && c0!=0) 
	    return (int)((-a[0]*a[4]+a[1]*a[3]-a[1]*y+a[4]*x)/c0);
	else if(a[2]==0.0 && a[1]!=0.0 && a[5]!=0.0)
	    return (int)((a[0]*a[4]-a[1]*a[3]+a[1]*y-a[4]*x)/(a[1]*a[5]));
	else
	    return -1;
*/
};
double  ezgdal_cr2y(EZGDAL_LAYER *layer, int col, int row) {
	double *a = ezgdal_layer_get_at(layer);
	return a[3]+a[4]*col+a[5]*row;
};




/*==========================================*/
/*                  LAYER                   */
/*                                          */

EZGDAL_LAYER*  ezgdal_open_layer(char *fname) {

  EZGDAL_LAYER *layer = NULL;
  
  if(!ezgdal_file_exists(fname)) {
    ezgdal_show_message(stderr,"Input file does not exist!");
    return NULL;
  }

  GDALAllRegister();

  layer = (EZGDAL_LAYER *)calloc(1,sizeof(EZGDAL_LAYER));

  layer->dataset_h = GDALOpen(fname,GA_ReadOnly);
  if(!(layer->dataset_h)) {
    ezgdal_show_message(stderr,"Problem with open input file!!");
    return NULL;
  }

  layer->band_h = GDALGetRasterBand(layer->dataset_h,1);

  layer->cols = GDALGetRasterBandXSize(layer->band_h);
  layer->rows = GDALGetRasterBandYSize(layer->band_h);

  layer->no_data = GDALGetRasterNoDataValue(layer->band_h, &(layer->is_no_data));
  layer->buffer = (double *)malloc(layer->cols*sizeof(double));
  layer->stripe = NULL;
  layer->frameset = NULL;

  return layer;
}

void  ezgdal_set_palette255(EZGDAL_LAYER *layer, double palette[][5], int n) {

  GDALColorTableH ct = GDALCreateColorTable(GPI_RGB);
  int i,i1,i2;
  GDALColorEntry c1,c2;

  for(i=0; i<n-1; i++) {
    c1.c1 = (short)palette[i][1];
    c1.c2 = (short)palette[i][2];
    c1.c3 = (short)palette[i][3];
    c1.c4 = (short)palette[i][4];
    i1    = (int)palette[i][0];
    c2.c1 = (short)palette[i+1][1];
    c2.c2 = (short)palette[i+1][2];
    c2.c3 = (short)palette[i+1][3];
    c2.c4 = (short)palette[i+1][4];
    i2    = (int)palette[i+1][0];
    GDALCreateColorRamp(ct,i1,&c1,i2,&c2);
  }
  if(palette[n-1][0]<255) {
    c1.c1 = 0;
    c1.c2 = 0;
    c1.c3 = 0;
    c1.c4 = 0;
    for(i=(int)palette[n-1][0]+1; i<=255; i++)
      GDALSetColorEntry(ct,i,&c1);
  }

  GDALSetRasterColorTable(layer->band_h,ct);
}


void  free_layer_stats(EZGDAL_LAYER *layer);


void  ezgdal_close_layer(EZGDAL_LAYER *layer) {

  GDALClose(layer->dataset_h);
  free(layer->buffer);
  free_layer_stats(layer);

  free(layer);
}

EZGDAL_LAYER*  ezgdal_create_layer(char *fname, 
                         char *wkt,
                         char *d_type,
                         double *geo_tr,
                         int rows,
                         int cols,
                         int *nodata) {

  EZGDAL_LAYER *o;
  int is_no_data = FALSE;
  double no_data = DBL_MIN;
  const char *data_format = "GTiff";
  GDALDataType data_type;
  GDALDriverH driver;
  
  o = (EZGDAL_LAYER *)calloc(1,sizeof(EZGDAL_LAYER));

  GDALAllRegister();

  driver = GDALGetDriverByName(data_format);
  if(driver==NULL) {
    ezgdal_show_message(stderr,"GDAL: Problem with GeoTiff driver!!");
    exit(1);
  }
    
  data_type = ezgdal_data_type(d_type);

  if(nodata) {
    no_data = *nodata;
    is_no_data = TRUE;
  }

  o->cols = cols;
  o->rows = rows;
  
  o->dataset_h = GDALCreate(driver,fname,o->cols,o->rows,1,data_type,NULL);

  if(!(o->dataset_h)) {
    ezgdal_show_message(stderr,"Problem with  file!!");
    exit(1);
  }
     
  o->band_h = GDALGetRasterBand(o->dataset_h,1);

  o->no_data = no_data;
  o->is_no_data = is_no_data;

  if(is_no_data)
    GDALSetRasterNoDataValue(o->band_h,no_data);

  o->buffer = (double *)malloc(o->cols*sizeof(double));
  o->stripe = NULL;
  o->frameset = NULL;

  GDALSetProjection(o->dataset_h,wkt);
  GDALSetGeoTransform(o->dataset_h,geo_tr);
  
  return o;
}

void  ezgdal_read_buffer(EZGDAL_LAYER *layer, int row) {
  int i;
  double d = 0.0;
  double *p;
  if(row<0 || row>=layer->rows) {
    if(layer->is_no_data) d = layer->no_data;
    p = layer->buffer;
    for(i=0; i<layer->cols; i++)
      *(p++) = d;
  } else
    i=GDALRasterIO(layer->band_h, GF_Read,0, row, layer->cols, 1,
                 layer->buffer, layer->cols, 1, GDT_Float64, 0, 0);
}

void  ezgdal_write_buffer(EZGDAL_LAYER *layer, int row) {
  if(row<0 || row>=layer->rows) return;
  CPLErr res = GDALRasterIO(layer->band_h, GF_Write, 0, row, layer->cols, 1,
               layer->buffer, layer->cols, 1, GDT_Float64, 0, 0);
  if(res>CE_Warning) {
    ezgdal_show_message(stderr,"GDAL I/O operation faild!");
    exit(EXIT_FAILURE);
  }
}

int  ezgdal_is_null(EZGDAL_LAYER *layer, double v) {
  if((layer->is_no_data) && (v == layer->no_data))
    return TRUE;
  return FALSE;
}

void  ezgdal_set_null(EZGDAL_LAYER *layer, double *v) {
  if(layer->is_no_data)
    *v = layer->no_data;
}



/*==========================================*/
/*                 STATS                    */
/*                                          */

void  ezgdal_calc_layer_stats(EZGDAL_LAYER *layer) {

  if(layer == NULL) return;

  if(layer->stats == NULL) {
    layer->stats = (EZGDAL_STATS *)malloc(sizeof(EZGDAL_STATS));

    assert(layer->stats!=NULL);

    layer->stats->map_cat = NULL;
    layer->stats->hist_N = 0;
    layer->stats->map_max_val = 0;
    layer->stats->hist_min = 0;
    layer->stats->hist_max = 0;
    layer->stats->hist_step = 0;
  }

  GDALComputeRasterStatistics(layer->band_h, FALSE,
                              &(layer->stats->min),
                              &(layer->stats->max),
                              &(layer->stats->avg),
                              &(layer->stats->std),
                              NULL,NULL);
}

void  ezgdal_calc_value_map(EZGDAL_LAYER *layer, double min, double max, int N) {

#ifdef OLD_GDAL
  int *hist;
#else
  GUIntBig *hist;
#endif

  int i,j = 0;
  CPLErr res;

  if(layer == NULL) return;
  if(min>=max || N<=0) return;

  if(layer->stats == NULL) {
    layer->stats = (EZGDAL_STATS *)malloc(sizeof(EZGDAL_STATS));

    assert(layer->stats!=NULL);

    layer->stats->min = 0;
    layer->stats->max = 0;
    layer->stats->avg = 0;
    layer->stats->std = 0;
  }

  layer->stats->hist_min = min;
  layer->stats->hist_max = max;
  layer->stats->hist_step = (double)N/(max-min);
  layer->stats->hist_N = N;

#ifdef OLD_GDAL
  hist = (int *) malloc(N * sizeof(int));
#else
  hist = (GUIntBig *) malloc(N * sizeof(GUIntBig));
#endif

#ifdef OLD_GDAL
  res = GDALGetRasterHistogram(layer->band_h,
                           min,
                           max,
                           N,
                           hist,
                           FALSE, FALSE,
                           NULL, NULL);
#else
  res = GDALGetRasterHistogramEx(layer->band_h,
                           min,
                           max,
                           N,
                           hist,
                           FALSE, FALSE,
                           NULL, NULL);
#endif

  assert(res==CE_None || res==CE_Debug);

  layer->stats->map_cat = (int *)malloc(N * sizeof(int));

  for(i=0; i<N; i++)
    if(hist[i]>0) 
      layer->stats->map_cat[i] = j++;
    else
      layer->stats->map_cat[i] = -1;

  layer->stats->map_max_val = j-1;

  free(hist);
}

int  ezgdal_get_value_index(EZGDAL_LAYER *layer, double val) {

  if(val < layer->stats->hist_min || val > layer->stats->hist_max)
    return -1;

  return layer->stats->map_cat[(int)( (val - layer->stats->hist_min) * layer->stats->hist_step )];
}

double  ezgdal_get_index_value(EZGDAL_LAYER *layer, int idx) {

  if(idx < 0 || idx >= layer->stats->hist_N)
    return -1;
    
  int i;
  for(i=0; i<layer->stats->hist_N; i++)
    if(layer->stats->map_cat[i]==idx)
      return (i + 0.5) / layer->stats->hist_step + layer->stats->hist_min;

  return -1;
}

void  free_layer_stats(EZGDAL_LAYER *layer) {

  if(layer == NULL) return;

  if(layer->stats != NULL) {
    if(layer->stats->map_cat != NULL)
      free(layer->stats->map_cat);
    free(layer->stats);
  }
}



/*==========================================*/
/*             STRIPE & FRAME               */
/*                                          */

EZGDAL_STRIPE*  ezgdal_create_stripe(EZGDAL_LAYER *layer, int row1, int height) {
  int r, c;
  EZGDAL_STRIPE *s;
  double d = 0.0, *p;

  if(layer==NULL) return NULL;

  if(row1<-height || row1>layer->rows-height-1) return NULL;

  if(layer->is_no_data) d = layer->no_data;
  s = malloc(sizeof(EZGDAL_STRIPE));
  s->layer = layer;
  s->rows = height;
  s->row1 = INT_MIN;
  s->row2 = s->row1 + height - 1;
  s->buffer = malloc((s->rows) * sizeof(double *));
  for(r=0; r<s->rows; r++) {
    s->buffer[r] = malloc((layer->cols + 2*height) * sizeof(double));
    p = s->buffer[r];
    for(c=0; c<layer->cols+2*height; c++)
      *(p++) = d;
  }
  s->frame = NULL;
  s->frames = 0;

  layer->stripe = s;

  return s;
}

void  reset_frame_buffer_pointers(EZGDAL_STRIPE *stripe) {
  int f,r;

  if(stripe->frame==NULL) return;

  for(f=0; f<stripe->frames; f++) {
    stripe->frame[f].row1 = stripe->row1;
    stripe->frame[f].row2 = stripe->row2;
    for(r=0; r<stripe->rows; r++) 
      stripe->frame[f].buffer[r] = stripe->buffer[r] + stripe->frame[f].col1;
  }
}

EZGDAL_FRAME*  ezgdal_create_frame(EZGDAL_STRIPE *stripe, int col1) {
  EZGDAL_FRAME *f;
  
  if(stripe == NULL) return NULL;
  if(stripe->frame!=NULL) return NULL;
  if(col1<-stripe->rows || col1>stripe->layer->cols-stripe->rows-1) return NULL;

  f = calloc(1,sizeof(EZGDAL_FRAME));
  f->owner.stripe = stripe;
  f->cols = stripe->rows;
  f->rows = stripe->rows;
  f->row1 = stripe->row1;
  f->row2 = stripe->row2;
  f->col1 = col1+stripe->rows;
  f->col2 = f->col1 + f->cols - 1;
  f->buffer = malloc(stripe->rows * sizeof(double *));
  stripe->frame = f;
  stripe->frames = 1;
  reset_frame_buffer_pointers(stripe);
  return f;
}

int  ezgdal_create_all_frames(EZGDAL_STRIPE *stripe, int start, int shift) {
  int r, s, N;

  if(stripe == NULL) return 0;
  if(stripe->frame!=NULL) return 0;

  if(start<-stripe->rows || start>stripe->layer->cols-stripe->rows-1) return 0;
  N = 0; 
  s = start + stripe->rows;
  while(s>=-stripe->rows && s<=stripe->layer->cols) {
    N++;
    s += shift;
  }
  stripe->frame = calloc(N, sizeof(EZGDAL_FRAME));
  stripe->frames = N;
  for(r=0; r<N; r++) {
    stripe->frame[r].owner.stripe = stripe;
    stripe->frame[r].cols = stripe->rows;
    stripe->frame[r].rows = stripe->rows;
    stripe->frame[r].row1 = stripe->row1;
    stripe->frame[r].row2 = stripe->row2;
    stripe->frame[r].col1 = start + stripe->rows;
    stripe->frame[r].col2 = stripe->frame[r].col1 + stripe->rows - 1;
    stripe->frame[r].buffer = malloc(stripe->rows * sizeof(double *));
    start += shift;
  }
  reset_frame_buffer_pointers(stripe);
  return N;
}

void  ezgdal_shift_frame_pos(EZGDAL_FRAME *frame, int col) {

  if(frame==NULL || frame->owner.stripe->frames!=1) return;
//  if(col<frame->cols || col>frame->owner.stripe->layer->cols-frame->cols-1) return;
  if(col<-frame->cols || col>frame->owner.stripe->layer->cols-frame->cols-1) return;

  frame->col1 = col + frame->cols;
  frame->col2 = frame->col1 + frame->cols - 1;
  reset_frame_buffer_pointers(frame->owner.stripe);
}

void  ezgdal_free_all_frames(EZGDAL_STRIPE *stripe) {
  int f;
  if(stripe->frame!=NULL) {
    for(f=0; f<stripe->frames; f++)
      free(stripe->frame[f].buffer);
    free(stripe->frame);
    stripe->frame = NULL;
    stripe->frames = 0;
  }
}

void  ezgdal_free_stripe(EZGDAL_LAYER *layer) {
  int r;
  ezgdal_free_all_frames(layer->stripe);
  for(r=0; r<layer->stripe->rows; r++)
    free(layer->stripe->buffer[r]);
  free(layer->stripe->buffer);
  free(layer->stripe);
  layer->stripe = NULL;
}

EZGDAL_FRAME*  ezgdal_get_frame(EZGDAL_STRIPE *stripe, int idx) {
  if(idx<0 || idx>=stripe->frames) return NULL;
  return &(stripe->frame[idx]);
}

int  ezgdal_load_stripe_data(EZGDAL_STRIPE *stripe, int row1) {
  int r, n, f;
  double *p;

  if(stripe == NULL) return 0;
  /* nothing to do */
  if(row1 == stripe->row1) return 0;

  /* out of range */
  if(row1<-stripe->rows || row1>stripe->layer->rows-1) return 0;

/*
    p=stripe->layer->buffer;
    for(r=0; r<stripe->rows; r++) {
      stripe->layer->buffer = stripe->buffer[r] + stripe->rows;
      ezgdal_read_buffer(stripe->layer, r + row1);
    }
    stripe->layer->buffer = p;
    n = stripe->rows;
return n;
*/

  if(row1 < stripe->row1 - stripe->rows ||
     row1 > stripe->row2) {
    /* read all rows */
    p=stripe->layer->buffer;
    for(r=0; r<stripe->rows; r++) {
      stripe->layer->buffer = stripe->buffer[r] + stripe->rows;
      ezgdal_read_buffer(stripe->layer,r + row1);
    }
    stripe->layer->buffer = p;
    n = stripe->rows;
  } else if(row1 > stripe->row1) {
    /* shift up and read only bottom */
    n = stripe->row2 -row1 + 1;
    for(r=0; r<n; r++) {
      p = stripe->buffer[r];
      stripe->buffer[r] = stripe->buffer[stripe->rows - n + r];
      stripe->buffer[stripe->rows - n + r] = p;
    }
    p=stripe->layer->buffer;
    for(r=n; r<stripe->rows; r++) {
      stripe->layer->buffer = stripe->buffer[r] + stripe->rows;
      ezgdal_read_buffer(stripe->layer,r + row1);
    }
    stripe->layer->buffer = p;
    reset_frame_buffer_pointers(stripe);
    n = stripe->rows-n;
  } else {
    /* shift down and read only top */
    n = row1 - stripe->row1 + 1;
    for(r=0; r<n; r++) {
      p = stripe->buffer[stripe->rows - n + r];
      stripe->buffer[stripe->rows - n + r] = stripe->buffer[r];
      stripe->buffer[r] = p;
    }
    p=stripe->layer->buffer;
    for(r=0; r<n; r++) {
      stripe->layer->buffer = stripe->buffer[r] + stripe->rows;
      ezgdal_read_buffer(stripe->layer,r + row1);
    }
    stripe->layer->buffer = p;
    reset_frame_buffer_pointers(stripe);
  }
  stripe->row1 = row1;
  stripe->row2 = row1 + stripe->rows - 1;

  for(f=0; f<stripe->frames; f++) {
    stripe->frame[f].row1 = stripe->row1;
    stripe->frame[f].row2 = stripe->row2;
  }

  return n;
}

int  ezgdal_save_stripe_data(EZGDAL_STRIPE *stripe) {
  int row1, N;
  double *p;

  if(stripe->row1<0) 
    row1 = 0;
  else
    row1 = stripe->row1;

  p=stripe->layer->buffer;
  N = 0;
  while(row1<=stripe->row2) {
    stripe->layer->buffer = stripe->buffer[row1] + stripe->rows;
    ezgdal_write_buffer(stripe->layer,row1);
    row1++;
    N++;
  }
  stripe->layer->buffer = p;
  return N;
}


/*==========================================*/
/*         FRAMESET & FRAME                 */
/*                                          */

void  ezgdal_unload_frameset_frame_data(EZGDAL_FRAME *frame) {
  if(frame->private_buffer!=NULL) {
    CPLFree(frame->private_buffer);
    free(frame->buffer);
    frame->private_buffer = NULL;
    frame->buffer = NULL;
  }
}

void  ezgdal_frameset_frame_alloc(EZGDAL_FRAME *frame) {
  int i;
  
  ezgdal_unload_frameset_frame_data(frame);

  unsigned long size = (unsigned long)(frame->col2-frame->col1+1)*(frame->row2-frame->row1+1);
  frame->private_buffer = (double *)CPLMalloc(size*sizeof(double));

  frame->buffer = (double **)malloc(frame->rows*sizeof(double *));

  if(frame->buffer==NULL || frame->private_buffer==NULL) {
    ezgdal_show_message(stderr,"No RAM to proceed!");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<frame->rows; i++) 
    frame->buffer[i] = frame->private_buffer + i*frame->cols;
}


EZGDAL_FRAMESET*  ezgdal_create_frameset_with_size(EZGDAL_LAYER *layer, int size) {
  EZGDAL_FRAMESET *fst = (EZGDAL_FRAMESET *)calloc(1,sizeof(EZGDAL_FRAMESET));
  fst->frame = (EZGDAL_FRAME **)calloc(size,sizeof(EZGDAL_FRAME *));
  fst->frameset_len = size;
  fst->frames = 0;
  fst->layer = layer;
  layer->frameset = fst;
  return fst;
}

EZGDAL_FRAMESET*  ezgdal_create_frameset(EZGDAL_LAYER *layer) {
  return ezgdal_create_frameset_with_size(layer,EZGDAL_FRAMESET_LEN);
}

void ezgdal_frameset_max_buffer_size(unsigned long size) {
  if(size>0) max_frame_buffer_size = size;
}

EZGDAL_FRAME*  ezgdal_add_frameset_frame(EZGDAL_FRAMESET *frameset, int col1, int col2, int row1, int row2) {
  int p;

  if(col1>col2) { p = col1; col1 = col2; col2 = p;}
  if(row1>row2) { p = col1; col1 = col2; col2 = p;}

  unsigned long size = (unsigned long)(col2-col1+1)*(row2-row1+1)*sizeof(double);

  if(size<=0 || size>max_frame_buffer_size)
    return NULL;

  EZGDAL_FRAME *frame = (EZGDAL_FRAME *)malloc(sizeof(EZGDAL_FRAME));
  frame->owner.frameset = frameset;
  frame->private_buffer = NULL;
  frame->buffer = NULL;

  ezgdal_frameset_set_frame(frame, col1, col2, row1, row2);
  
  if(frameset->frames==frameset->frameset_len) {
    frameset->frameset_len += EZGDAL_FRAMESET_STEP;
    frameset->frame = realloc(frameset->frame, frameset->frameset_len);
    if(frameset->frame==NULL) {
      ezgdal_show_message(stderr,"No RAM to proceed!");
      exit(EXIT_FAILURE);
    }
  }
  frameset->frame[frameset->frames] = frame;
  frameset->frames++;
  
  return frame;
}

void  ezgdal_frameset_set_frame(EZGDAL_FRAME *frame, int col1, int col2, int row1, int row2) {
  int p;
  int data_are_loaded;
  
  if(frame==NULL) return;

  data_are_loaded = (frame->private_buffer!=NULL);
  ezgdal_unload_frameset_frame_data(frame);

  if(col1>col2) { p = col1; col1 = col2; col2 = p;}
  if(row1>row2) { p = row1; row1 = row2; row2 = p;}

  frame->col1 = col1;
  frame->col2 = col2;
  frame->row1 = row1;
  frame->row2 = row2;
  frame->cols = col2-col1+1;
  frame->rows = row2-row1+1;
  
  if(data_are_loaded)
    ezgdal_load_frameset_frame_data(frame);
}

void  ezgdal_free_frameset(EZGDAL_FRAMESET *frameset) {
  if(frameset==NULL) return;
  EZGDAL_LAYER *l = frameset->layer;
  ezgdal_free_frameset_all_frames(frameset);
  free(frameset->frame);
  free(frameset);
  l->frameset = NULL;
}

void  ezgdal_free_frameset_frame(EZGDAL_FRAME *frame) {
  ezgdal_unload_frameset_frame_data(frame);
  free(frame);
}

void  ezgdal_free_frameset_all_frames(EZGDAL_FRAMESET *frameset) {
  int i;
  for(i=0; i<frameset->frames; i++)
    ezgdal_free_frameset_frame(frameset->frame[i]);
  frameset->frames = 0;
}

EZGDAL_FRAME*  ezgdal_get_frameset_frame(EZGDAL_FRAMESET *frameset, int idx) {
  EZGDAL_FRAME *frame = NULL;
  if(idx<frameset->frames && idx>=0)
    frame = frameset->frame[idx];
  return frame;
}

void  ezgdal_load_frameset_frame_data(EZGDAL_FRAME *frame) {

  if(frame==NULL) return;
  ezgdal_frameset_frame_alloc(frame);

  EZGDAL_LAYER *l = frame->owner.frameset->layer;

  if(frame->col1 >= 0 &&
     frame->row1 >= 0 &&
     frame->col2 < l->cols &&
     frame->row2 < l->rows) {
    // inside - read
    CPLErr res = GDALRasterIO(l->band_h,
                                GF_Read, 
                                frame->col1, frame->row1,
                                frame->cols, frame->rows,
                                frame->private_buffer, 
                                frame->cols, frame->rows,
                                GDT_Float64, 0, 0 );
      
    if(res>CE_Warning) {
      ezgdal_show_message(stderr,"GDAL I/O operation faild!");
      exit(EXIT_FAILURE);
    }
    
  } else {
    if(!(frame->col1 >= l->cols ||
         frame->row1 >= l->rows ||
         frame->col2 < 0 || frame->row2 < 0)) {
      // border crossing - read part and copy

      int row, col, r, c, new_c1, new_r1, new_c2, new_r2, new_cols, new_rows;

      int i;
      long n;
      double *p;

      n = frame->cols*frame->rows;
      double d = 0.0;
      if(l->is_no_data)
        d = l->no_data;
      else if(l->stats!=NULL) {
        if(l->stats->min<=0.0 && l->stats->max>=0.0)
          d = 0.0;
        else 
          d = l->stats->min;
      }
      p = frame->private_buffer;
      for(i=0; i<n; i++)
        *(p++) = d;

      new_c1 = (frame->col1<0)?0:frame->col1;
      new_r1 = (frame->row1<0)?0:frame->row1;
      new_c2 = (frame->col2>=l->cols)?l->cols-1:frame->col2;
      new_r2 = (frame->row2>=l->rows)?l->rows-1:frame->row2;
      new_cols = new_c2 - new_c1 + 1;
      new_rows = new_r2 - new_r1 + 1;
      double *buf = (double *)malloc(new_cols*new_rows*sizeof(double));
      
      CPLErr res = GDALRasterIO(l->band_h,
                                GF_Read, 
                                new_c1, new_r1,
                                new_cols, new_rows,
                                buf, 
                                new_cols, new_rows,
                                GDT_Float64, 0, 0);
      
      if(res>CE_Warning) {
        free(buf);
        ezgdal_show_message(stderr,"GDAL I/O operation faild!");
        exit(EXIT_FAILURE);
      }

      for(r=0; r<new_rows; r++) {
        row = new_r1 - frame->row1 + r;
        col = new_c1 - frame->col1;
        for(c=0; c<new_cols; c++) 
          frame->private_buffer[row*frame->cols+col+c] = buf[r*new_cols+c];
        row += frame->cols;
      }

      free(buf);
    }
  }

}

