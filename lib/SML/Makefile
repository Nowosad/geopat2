
PROG = libsml

COMMON =  smllib.c smltxt.c
COMMON_O = $(COMMON:%.c=%.o)
HEADERS = sml.h

all install: $(PROG).h $(PROG).a

$(PROG).a: $(COMMON_O)
	$(AR) $(PROG).a $(COMMON_O)

$(PROG).h: $(HEADERS)
	cat $(HEADERS) > $(PROG).h

$(COMMON_O): $(COMMON) $(HEADERS)
	$(CC) $(CFLAGS) $(LIBCFLAGS) -c  $(COMMON)

clean:
	rm -f *.o $(PROG).h $(PROG).a
