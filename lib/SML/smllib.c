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


#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#define DLL_EXPORT

#include "sml.h"

#define _FILE_OFFSET_BITS 64

void error(char *err) {
	if(err!=NULL) fprintf(stderr,err);
	if(errno!=0)
		fprintf(stderr,"\nSML error: %s\n",strerror(errno));
	fflush(stderr);
	exit(0);
}


long int sml_get_file_size(FILE *f) {
    long int p,i;
    p = ftell(f);
    fseek(f,0,SEEK_END);
    i = ftell(f);
    fseek(f,p,SEEK_SET);
    return i;
}

SML_WINDOW *sml_create_window(int rows, int cols, 
		    double at0, double at1, double at2, 
		    double at3, double at4, double at5, 
		    char* proj) {
	SML_WINDOW *win = (SML_WINDOW *)calloc(1,sizeof(SML_WINDOW));
	win->at[0] = at0;
	win->at[1] = at1;
	win->at[2] = at2;
	win->at[3] = at3;
	win->at[4] = at4;
	win->at[5] = at5;
	win->rows = rows;
	win->cols = cols;
	strncpy(win->proj,proj,SML_DATA_PROJ_LEN-1);
	return win;
}

SML_WINDOW *sml_create_window_copy(SML_WINDOW *w) {
	int i;
	SML_WINDOW *win = (SML_WINDOW *)calloc(1,sizeof(SML_WINDOW));
	for(i=0; i<6; i++)
	  win->at[i] = w->at[i];
	win->rows = w->rows;
	win->cols = w->cols;
	strncpy(win->proj,w->proj,SML_DATA_PROJ_LEN-1);
	return win;
}

void sml_free_window(SML_WINDOW *w) {
	free(w);
}

SML_CELL_TYPE *sml_create_cell_type(SML_D_TYPE d_type, int ndim, int* dims) {
	int i;
	SML_CELL_TYPE *ct = (SML_CELL_TYPE *)malloc(sizeof(SML_CELL_TYPE));
	ct->d_type = d_type;
	ct->dim = ndim;
	if(ndim>0) {
		ct->dims = (int *)malloc(sizeof(int)*ndim);
		ct->len = (int *)malloc(sizeof(int)*ndim);
		for(i=0; i<ndim; i++) 
		  ct->dims[i] = dims[i];
		ct->len[ndim-1] = 1;
		for(i=ndim-2; i>=0; i--)
			ct->len[i]=ct->len[i+1]*ct->dims[i];
	} else {
		ct->dims = NULL;
		ct->len = NULL;
	}
	return ct;
}

SML_CELL_TYPE *sml_create_cell_type_copy(SML_CELL_TYPE *ct0) {
	int i;
	SML_CELL_TYPE *ct = (SML_CELL_TYPE *)malloc(sizeof(SML_CELL_TYPE));
	ct->d_type = ct0->d_type;
	ct->dim = ct0->dim;
	if(ct0->dim>0) {
		ct->dims = (int *)malloc(sizeof(int)*ct0->dim);
		ct->len = (int *)malloc(sizeof(int)*ct0->dim);
		for(i=0; i<ct0->dim; i++) 
		  ct->dims[i] = ct0->dims[i];
		ct->len[ct0->dim-1] = 1;
		for(i=ct0->dim-2; i>=0; i--)
			ct->len[i]=ct->len[i+1]*ct->dims[i];
	} else {
		ct->dims = NULL;
		ct->len = NULL;
	}
	return ct;
}

void sml_free_cell_type(SML_CELL_TYPE *ct) {
	if(ct->dim > 0) {
		free(ct->dims);
		free(ct->len);
	}
	free(ct);
}

void sml_recalc_data_header(SML_DATA_HEADER *dh) {
	int i,j,k;

	j=1;
	if(dh->cell_type->dim > 0) 
		for(i=0; i<dh->cell_type->dim; i++) 
			j*=dh->cell_type->dims[i];
	dh->cell_N_elements = j;
	dh->cell_type->len=(int *)malloc(sizeof(int)*dh->cell_type->dim);
	dh->cell_type->len[dh->cell_type->dim-1] = 1;
		for(i=dh->cell_type->dim-2; i>=0; i--)
		dh->cell_type->len[i]=dh->cell_type->len[i+1]*dh->cell_type->dims[i];
	k=0;
	switch(dh->cell_type->d_type) {
		case SML_BYTE: 
				k=sizeof(char);
				break;
		case SML_SHORT_INT: 
				k=sizeof(short int);
				break;
		case SML_INT: 
				k=sizeof(int);
				break;
		case SML_FLOAT: 
				k=sizeof(float);
				break;
		case SML_DOUBLE:
				k=sizeof(double);
				break;
	}
	dh->cell_element_size=k;
	dh->cell_size=j*k+1;
}

int sml_write_layer_header(SML_DATA_HEADER *dh) {
	char name[512];
	FILE *f;
	int i;
	
	strcpy(name,dh->name);
	strcat(name,".hdr");
	f=fopen(name,"w");
	fprintf(f,"dim: %d\n",dh->cell_type->dim);
	if(dh->cell_type->dim > 0) {
		fprintf(f,"dims: %d",dh->cell_type->dims[0]);
		for(i=1; i < dh->cell_type->dim; i++)
			fprintf(f,",%d",dh->cell_type->dims[i]);
		fprintf(f,"\n");
	}
	switch(dh->cell_type->d_type) {
		case SML_DOUBLE:
				fprintf(f,"type: DOUBLE\n");
				break;
		case SML_FLOAT:
				fprintf(f,"type: FLOAT\n");
				break;
		case SML_INT:
				fprintf(f,"type: INT\n");
				break;
		case SML_SHORT_INT:
				fprintf(f,"type: SHORT_INT\n");
				break;
		case SML_BYTE:
				fprintf(f,"type: BYTE\n");
				break;
	}
	fprintf(f,"at0: %.18lf\n",dh->file_win->at[0]);
	fprintf(f,"at1: %.18lf\n",dh->file_win->at[1]);
	fprintf(f,"at2: %.18lf\n",dh->file_win->at[2]);
	fprintf(f,"at3: %.18lf\n",dh->file_win->at[3]);
	fprintf(f,"at4: %.18lf\n",dh->file_win->at[4]);
	fprintf(f,"at5: %.18lf\n",dh->file_win->at[5]);
	fprintf(f,"rows: %d\n",dh->file_win->rows);
	fprintf(f,"cols: %d\n",dh->file_win->cols);
	fprintf(f,"proj: %s\n",dh->file_win->proj);
	fprintf(f,"desc: %s\n",dh->desc);
	fclose(f);

	return 0;
}

int sml_read_layer_header(SML_DATA_HEADER *dh) {
	char name[4096];
	FILE *f;
	int i;
	strcpy(name,dh->name);
	strcat(name,".hdr");
	f=fopen(name,"r");
	if(!f) error("\nHeader file does not exist\n");
	
	dh->cell_type=(SML_CELL_TYPE *)malloc(sizeof(SML_CELL_TYPE));
	if(fscanf(f,"%*[^:]: %d\n",&(dh->cell_type->dim))!=1) error("\nError in reading header file\n");
	if(dh->cell_type->dim > 0) {
		dh->cell_type->dims=(int *)malloc(sizeof(int)*dh->cell_type->dim);
		if(fscanf(f,"%*[^:]: %d",&(dh->cell_type->dims[0]))!=1) error("\nError in reading header file\n");
		for(i=1; i < dh->cell_type->dim; i++)
			if(fscanf(f,",%d",&(dh->cell_type->dims[i]))!=1) error("\nError in reading header file\n");
	}
	if(fscanf(f,"%*[^: ]: %s\n",name)!=1) error("\nError in reading header file\n");
	if(strcmp(name,"DOUBLE")==0)
		dh->cell_type->d_type=SML_DOUBLE;
	else if(strcmp(name,"FLOAT")==0)
		dh->cell_type->d_type=SML_FLOAT;
	else if(strcmp(name,"INT")==0)
		dh->cell_type->d_type=SML_INT;
	else if(strcmp(name,"SHORT_INT")==0)
		dh->cell_type->d_type=SML_SHORT_INT;
	else if(strcmp(name,"BYTE")==0)
		dh->cell_type->d_type=SML_BYTE;
	else
		error("\nError in reading header file\n");
	dh->file_win = (SML_WINDOW *)calloc(1,sizeof(SML_WINDOW));
	if(fscanf(f,"%*[^:]: %lf",&(dh->file_win->at[0]))!=1) error("\nError in reading header file\n");
	if(fscanf(f,"%*[^:]: %lf",&(dh->file_win->at[1]))!=1) error("\nError in reading header file\n");
	if(fscanf(f,"%*[^:]: %lf",&(dh->file_win->at[2]))!=1) error("\nError in reading header file\n");
	if(fscanf(f,"%*[^:]: %lf",&(dh->file_win->at[3]))!=1) error("\nError in reading header file\n");
	if(fscanf(f,"%*[^:]: %lf",&(dh->file_win->at[4]))!=1) error("\nError in reading header file\n");
	if(fscanf(f,"%*[^:]: %lf",&(dh->file_win->at[5]))!=1) error("\nError in reading header file\n");
	if(fscanf(f,"%*[^:]: %d",&(dh->file_win->rows))!=1) error("\nError in reading header file\n");
	if(fscanf(f,"%*[^:]: %d",&(dh->file_win->cols))!=1) error("\nError in reading header file\n");

	fscanf(f,"%*[^:]: %[^\n]",dh->file_win->proj);
	fscanf(f,"%*[^:]: %[^\n]",dh->desc);
/*
	if(fscanf(f,"%*[^:]: %[^\n]",dh->file_win->proj)!=1) error("\nError in reading header file\n");
	
	if(fscanf(f,"%*[^:]: %[^\n]",dh->desc)!=1) error("\nError in reading header file\n");
*/
	fclose(f);

	return 0;
}

SML_DATA_HEADER* sml_create_layer(
			char* fname, 
			SML_CELL_TYPE *cell_type,
			SML_WINDOW *w) {

	SML_DATA_HEADER *dh = (SML_DATA_HEADER *)calloc(1,sizeof(SML_DATA_HEADER));
	dh->f = fopen(fname,"wb");
	strcpy(dh->name,fname);
	dh->file_status=SML_NEW;

	dh->cell_type = cell_type;
	dh->file_win = w;

	sml_recalc_data_header(dh);

	return dh;
};

SML_DATA_HEADER* sml_open_layer(char* fname) {

	SML_DATA_HEADER *dh = (SML_DATA_HEADER *)calloc(1,sizeof(SML_DATA_HEADER));
	dh->f = fopen(fname,"rb");
	if(dh->f == NULL) error("\nData file can not be opened\n");
	strcpy(dh->name,fname);
	dh->file_status = SML_EXISTING;

	sml_read_layer_header(dh);

	sml_recalc_data_header(dh);

	return dh;
};

void sml_close_layer(SML_DATA_HEADER *dh) {
	fclose(dh->f);
	if(dh->file_status == SML_NEW)
	  sml_write_layer_header(dh);
	sml_free_cell_type(dh->cell_type);
	sml_free_window(dh->file_win);
	free(dh);
};

void sml_set_layer_description(SML_DATA_HEADER *dh, char *desc[], int cnt) {
    char *s;
    int i,l;

    l = SML_DATA_DESC_LEN-1;
    s = dh->desc;
    *s = 0;
    for(i = 0; i < cnt; i++) {
	strncpy(s,desc[i],l);
	l-=(int)strlen(desc[i]);
        if(l<1) break;
        s+=strlen(desc[i]);
        if(i<cnt-1) {
	    strncpy(s,"|",l--);
    	    if(l<1) break;
    	    s+=1;
        }
    }
    *s=0;
}


void *sml_create_cell_row_buffer(SML_DATA_HEADER *dh) {
  unsigned int size = dh->file_win->cols * dh->cell_size;
  return malloc(size);
}

void *sml_get_cell_pointer(SML_DATA_HEADER *dh, void *row_buffer, int i) {
  void *p = NULL;
  if(i<dh->file_win->cols)
    p = (char *)row_buffer + i * dh->cell_size;
  return p;
}


void sml_read_cell_from_layer(SML_DATA_HEADER *dh, void *cell, int col, int row) {
    if(col<0 || col>=dh->file_win->cols || row<0 || row>=dh->file_win->rows) {
      sml_set_cell_null(cell);
      return;
    }
    unsigned long pos = ((unsigned long)row*(unsigned long)(dh->file_win->cols)+(unsigned long)col)*(unsigned long)(dh->cell_size);
    if(fseek(dh->f,pos,SEEK_SET)!=0) error("");
    if(fread(cell,dh->cell_size,1,dh->f)!=1) error("");
}

void sml_read_cell_from_layer_xy(SML_DATA_HEADER *dh, void *cell, double x, double y) {
    int row = sml_xy2r(dh,x,y);
    int col = sml_xy2c(dh,x,y);
    sml_read_cell_from_layer(dh,cell,col,row);
}

void sml_read_row_from_layer(SML_DATA_HEADER *dh, void *cell_row, int row) {
    unsigned long pos = (unsigned long)row*(unsigned long)(dh->file_win->cols)*(unsigned long)(dh->cell_size);
    if(fseek(dh->f,pos,SEEK_SET)!=0) error("");
    if(fread(cell_row,dh->cell_size,dh->file_win->cols,dh->f)!=dh->file_win->cols)
        error("");
}

void sml_write_row_to_layer(SML_DATA_HEADER *dh, void *cell_row, int row) {
    unsigned long pos = (unsigned long)row*(unsigned long)(dh->file_win->cols)*(unsigned long)(dh->cell_size);
    if(fseek(dh->f,pos,SEEK_SET)!=0) error("");
    if(fwrite(cell_row,dh->cell_size,dh->file_win->cols,dh->f)!=dh->file_win->cols)
        error("");
}

void sml_write_next_row_to_layer(SML_DATA_HEADER *dh, void *cell_row) {
    if(fwrite(cell_row,dh->cell_size,dh->file_win->cols,dh->f)!=dh->file_win->cols)
        error("");
}



int sml_xy2c(SML_DATA_HEADER *dh, double x, double y) {
	double *a = dh->file_win->at;
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
double sml_cr2x(SML_DATA_HEADER *dh, int col, int row) { 
	double *a = dh->file_win->at;
	return a[0]+a[1]*col+a[2]*row;
};

int sml_xy2r(SML_DATA_HEADER *dh, double x, double y) {
	double *a = dh->file_win->at;
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
double sml_cr2y(SML_DATA_HEADER *dh, int col, int row) {
	double *a = dh->file_win->at;
	return a[3]+a[4]*col+a[5]*row;
};


char sml_is_cell_null(void *cell) {
	assert(cell!=NULL);
	return *((char *)cell)==0;
};
void sml_set_cell_null(void *cell) {
	assert(cell!=NULL);
	*((char *)cell) = 0;
};
void sml_set_cell_not_null(void *cell) {
	assert(cell!=NULL);
	*((char *)cell) = 1;
};


void *sml_get_cell_data(void *cell) {
	return (char *)cell+1;
};

double sml_get_cell_val_dbl(SML_DATA_HEADER *dh, void *cell, int i) {
	double v = 0.0;
	if(i<dh->cell_N_elements) {
		int pos = i*dh->cell_element_size;
		switch(dh->cell_type->d_type) {
			case SML_DOUBLE:
				v=*((double *)((char *)sml_get_cell_data(cell)+pos));
				break;
			case SML_FLOAT:
				v=*((float *)((char *)sml_get_cell_data(cell)+pos));
				break;
			case SML_INT:
				v=*((int *)((char *)sml_get_cell_data(cell)+pos));
				break;
			case SML_SHORT_INT:
				v=*((short int *)((char *)sml_get_cell_data(cell)+pos));
				break;
			case SML_BYTE:
				v=*((char *)((char *)sml_get_cell_data(cell)+pos));
				break;
		}
	}
	return v;
};
int sml_get_cell_val_int(SML_DATA_HEADER *dh, void *cell, int i) {
	int v = 0;
	int pos = i*dh->cell_element_size;
	if(pos<dh->cell_size) {
		switch(dh->cell_type->d_type) {
			case SML_DOUBLE:
				v=(int)*((double *)((char *)sml_get_cell_data(cell)+pos));
				break;
			case SML_FLOAT:
				v=(int)*((float *)((char *)sml_get_cell_data(cell)+pos));
				break;
			case SML_INT:
				v=*((int *)((char *)sml_get_cell_data(cell)+pos));
				break;
			case SML_SHORT_INT:
				v=(int)*((short int *)((char *)sml_get_cell_data(cell)+pos));
				break;
			case SML_BYTE:
				v=(int)*((char *)((char *)sml_get_cell_data(cell)+pos));
				break;
		}
	}
	return v;
};
void *sml_get_cell_val(SML_DATA_HEADER *dh, void *cell, int i) {
	void *v = NULL;
	int pos = i*dh->cell_element_size;
	if(pos<dh->cell_size) 
		v=(char *)sml_get_cell_data(cell)+pos;
	return v;
};


void sml_set_cell_val_dbl(SML_DATA_HEADER *dh, double v, void *cell, int i) {
	if(i<dh->cell_N_elements) {
		int pos = i*dh->cell_element_size;
		int vv;
		switch(dh->cell_type->d_type) {
			case SML_DOUBLE:
				*((double *)((char *)sml_get_cell_data(cell)+pos))=v;
				break;
			case SML_FLOAT:
				*((float *)((char *)sml_get_cell_data(cell)+pos))=(float)v;
				break;
			case SML_INT:
				vv=(int)v;
				if(vv<-2147483647) vv=-2147483647;
				else if(vv>2147483647) vv=2147483647;
				*((int *)((char *)sml_get_cell_data(cell)+pos))=vv;
				break;
			case SML_SHORT_INT:
				vv=(int)v;
				if(vv<-32767) vv=-32767;
				else if(vv>32767) vv=-32767;
				*((short int *)((char *)sml_get_cell_data(cell)+pos))=vv;
				break;
			case SML_BYTE:
				vv=(int)v;
				if(vv<0) vv=0;
				else if(vv>255) vv=255;
				*((unsigned char *)((char *)sml_get_cell_data(cell)+pos))=vv;
				break;
		}
	}
};
void sml_set_cell_val_int(SML_DATA_HEADER *dh, int v, void *cell, int i) {
	int pos = i*dh->cell_element_size;
	if(pos<dh->cell_size) {
		switch(dh->cell_type->d_type) {
			case SML_DOUBLE:
				*((double *)((char *)sml_get_cell_data(cell)+pos))=v;
				break;
			case SML_FLOAT:
				*((float *)((char *)sml_get_cell_data(cell)+pos))=(float)v;
				break;
			case SML_INT:
				*((int *)((char *)sml_get_cell_data(cell)+pos))=v;
				break;
			case SML_SHORT_INT:
				if(v<-32767) v=-32767;
				else if(v>32767) v=-32767;
				*((short int *)((char *)sml_get_cell_data(cell)+pos))=v;
				break;
			case SML_BYTE:
				if(v<0) v=0;
				else if(v>255) v=255;
				*((char *)((char *)sml_get_cell_data(cell)+pos))=v;
				break;
		}
	}
};
