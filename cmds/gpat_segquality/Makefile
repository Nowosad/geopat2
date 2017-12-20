
PROG = gpat_segquality

LIBS = \
../../lib/measures/libmeasures.a \
../../lib/tools/libtools.a \
../../lib/argtable/libargtable3.a

FILES = \
compatibility.c \
list.c \
hexagonal.c \
calculate.c \
write_results.c \
main.c \
diagnostics.c \
rebuild.c

OFILES = $(FILES:%.c=%.o)

ifndef PREFIX
	PREFIX = /usr/local
endif

all: $(PROG)

%.o : %.c
	$(CC) $(CFLAGS) $(EXTFLAGS) -c $<

$(PROG): $(OFILES)
	$(CC) -o $(PROG) $(OFILES) $(LIBS) $(CFLAGS) $(EXTFLAGS)

install: $(PROG)
	mv -f $(PROG) $(PREFIX)/bin

clean:
	rm -f $(PROG) *.o
