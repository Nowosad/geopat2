CC = gcc
CFLAGS = -Wall -I/usr/include/gdal -fopenmp -O2 -I/usr/local/include -DSML_LINUX -DEZGDAL_LINUX
LIBCFLAGS = 
EXTFLAGS = -lm -lgomp -lgdal -L/usr/local/lib -lsml -lezgdal 
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
