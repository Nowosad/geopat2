
PROG = gpat_segment

LIBS = \
../../lib/signatures/libsignatures.a \
../../lib/measures/libmeasures.a \
../../lib/tools/libtools.a \
../../lib/argtable/libargtable3.a

FILES = \
compatibility.c \
list.c \
cluster_growing.c \
cluster_heterogeneity.c \
cluster_hierarhical.c \
cluster_minarea.c \
cluster_seeds.c \
hexagonal.c \
calculate.c \
post_swap.c \
write_results.c \
main.c

OFILES = $(FILES:%.c=%.o)

ifndef PREFIX
	PREFIX = /usr/local
endif

all: $(PROG)

%.o : %.c
	$(CC) $(CFLAGS) $(EXTFLAGS) -c $<

$(PROG): $(OFILES)
	$(CC) $(CFLAGS) $(EXTFLAGS) -o $(PROG) $(OFILES) $(LIBS)

install: $(PROG)
	mv -f $(PROG) $(PREFIX)/bin

clean:
	rm -f $(PROG) *.o
