CC = gcc
CFLAGS = -Wall -I/usr/include/gdal -fopenmp -O2 -DSML_LINUX -DEZGDAL_LINUX
LIBCFLAGS = 
EXTFLAGS = -lm -lgomp -L/usr/lib/gdal -lgdal -lsml -lezgdal 
AR = ar cvq

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
