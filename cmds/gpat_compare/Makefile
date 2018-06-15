
PROG = gpat_compare

LIBS = \
../../lib/SML/libsml.a \
../../lib/ezGDAL/libezgdal.a \
../../lib/signatures/libsignatures.a \
../../lib/measures/libmeasures.a \
../../lib/tools/libtools.a \
../../lib/argtable/libargtable3.a

ifndef PREFIX
	PREFIX = /usr/local
endif



all: $(PROG)

$(PROG): main.o palette.o
	$(CC) -o $(PROG)  main.o palette.o $(LIBS) $(CFLAGS) $(EXTFLAGS) 

main.o: main.c ../../lib/argtable/argtable3.h palette.h
	$(CC) $(CFLAGS) $(EXTFLAGS) -c main.c

palette.o: palette.c palette.h
	$(CC) $(CFLAGS) $(EXTFLAGS) -c palette.c

install: $(PROG)
	mv -f $(PROG) $(PREFIX)/bin

clean:
	rm -f $(PROG) *.o
