
SUBDIRS= \
argtable \
measures \
normalization \
signatures \
tools

all clean install:
	@for dir in $(SUBDIRS); do \
	    ( \
		echo; \
		echo ================================; \
		echo   Library: $$dir; \
		echo ================================; \
		cd $$dir; \
		$(MAKE) $@; \
		); \
	done
