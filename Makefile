CC = gcc
CFLAGS = -Wall -g -I/usr/include/gdal -fopenmp -O2 -DSML_LINUX -DEZGDAL_LINUX
LIBCFLAGS = 
EXTFLAGS = -lm -lgomp -lgdal
AR = ar cvq

ifndef PREFIX
  PREFIX = /usr/local
endif

export CC
export AR
export CFLAGS
export LIBCFLAGS
export EXTFLAGS
export PREFIX

SUBDIRS= \
lib \
cmds 

# wgui

all clean install:
	@for dir in $(SUBDIRS); do \
	    (cd $$dir; $(MAKE) $@); \
	done
