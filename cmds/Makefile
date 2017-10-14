
SUBDIRS= \
gpat_compare \
gpat_distmtx \
gpat_globnorm \
gpat_grd2txt \
gpat_gridhis \
gpat_gridts \
gpat_pointshis \
gpat_pointsts \
gpat_polygon \
gpat_search \
gpat_segment \
gpat_segquality

all clean install:
	@for dir in $(SUBDIRS); do \
	    ( \
		echo; \
		echo ================================; \
		echo   $$dir; \
		echo ================================; \
		cd $$dir; \
		$(MAKE) $@; \
		); \
	done
