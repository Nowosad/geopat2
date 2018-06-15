
PROG = libmeasures

MEASURES = $(shell ls measure_*.c)
MEASURES_O = $(MEASURES:%.c=%.o)
COMMON = measures_interface.c geopat_compatibility.c emd.c
COMMON_O = $(COMMON:%.c=%.o)
HEADERS = $(shell ls *.h)

all install: $(PROG).a

$(PROG).a: $(MEASURES_O) $(COMMON_O)
	$(AR) $(PROG).a $(MEASURES_O) $(COMMON_O)

$(MEASURES_O): $(MEASURES) $(HEADERS)
	$(CC) $(CFLAGS) -c $(MEASURES)

measures_interface.o: measures_interface.c
	$(CC) $(CFLAGS) $(LIBCFLAGS) -c measures_interface.c
	
geopat_compatibility.o: geopat_compatibility.c
	$(CC) $(CFLAGS) $(LIBCFLAGS) -c geopat_compatibility.c
	
emd.o: emd.c
	$(CC) $(CFLAGS) $(LIBCFLAGS) -c emd.c

clean:
	rm -f *.o $(PROG).a
