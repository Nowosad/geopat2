
PROG = libsignatures

SIGNATURES = $(shell ls signature_*.c)
SIGNATURES_O = $(SIGNATURES:%.c=%.o)
COMMON = signatures_interface.c
COMMON_O = $(COMMON:%.c=%.o)
HEADERS = $(shell ls *.h)


all install: $(PROG).a

$(PROG).a: $(SIGNATURES_O) $(COMMON_O)
	$(AR) $(PROG).a $(SIGNATURES_O) $(COMMON_O)

$(SIGNATURES_O): $(SIGNATURES) $(HEADERS)
	$(CC) -c $(CFLAGS) $(LIBCFLAGS) $(SIGNATURES)

$(COMMON_O): $(COMMON) $(HEADERS)
	$(CC) -c $(CFLAGS) $(LIBCFLAGS) $(COMMON)

clean:
	rm -f *.o $(PROG).a
