
PROG = gpat_grd2txt

LIBS = \
../../lib/tools/libtools.a \
../../lib/argtable/libargtable3.a

ifndef PREFIX
	PREFIX = /usr/local
endif

all: $(PROG)

$(PROG): main.o
	$(CC) $(CFLAGS) $(EXTFLAGS) -o $(PROG) main.o $(LIBS)

main.o: main.c ../../lib/argtable/argtable3.h
	$(CC) $(CFLAGS) $(EXTFLAGS) -c main.c

install: $(PROG)
	mv -f $(PROG) $(PREFIX)/bin

clean:
	rm -f $(PROG) *.o