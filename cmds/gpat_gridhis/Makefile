
PROG = gpat_gridhis

LIBS = \
../../lib/SML/libsml.a \
../../lib/ezGDAL/libezgdal.a \
../../lib/signatures/libsignatures.a \
../../lib/normalization/libnormalization.a \
../../lib/tools/libtools.a \
../../lib/argtable/libargtable3.a

ifndef PREFIX
	PREFIX = /usr/local
endif

all: $(PROG)

$(PROG): main.o
	$(CC) -o $(PROG) main.o $(LIBS) $(CFLAGS) $(EXTFLAGS)

main.o: main.c
	$(CC) $(CFLAGS) $(EXTFLAGS) -c main.c

install: $(PROG)
	mv -f $(PROG) $(PREFIX)/bin

clean:
	rm -f $(PROG) *.o
