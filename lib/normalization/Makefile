
PROG = libnormalization

METHODS = $(shell ls method_*.c)
METHODS_O = $(METHODS:%.c=%.o)
COMMON = methods_interface.c
COMMON_O = $(COMMON:%.c=%.o)
HEADERS = $(shell ls *.h)


all install: $(PROG).a

$(PROG).a: $(METHODS_O) $(COMMON_O)
	$(AR) $(PROG).a $(METHODS_O) $(COMMON_O)

$(METHODS_O): $(METHODS) $(HEADERS)
	$(CC) -c $(CFLAGS) $(LIBCFLAGS) $(METHODS)

$(COMMON_O): $(COMMON) $(HEADERS)
	$(CC) -c $(CFLAGS) $(LIBCFLAGS) $(COMMON)

clean:
	rm -f *.o $(PROG).a
