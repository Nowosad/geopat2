
PROG = libargtable3

COMMON = argtable3.c
COMMON_O = $(COMMON:%.c=%.o)
HEADERS = argtable3.h

all install: $(PROG).a

$(PROG).a: $(COMMON_O)
	$(AR) $(PROG).a $(COMMON_O)

$(COMMON_O): $(COMMON) $(HEADERS)
	$(CC) $(CFLAGS) $(LIBCFLAGS) -c $(COMMON)

clean:
	rm -f *.o $(PROG).a
