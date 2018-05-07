
PROG = libtools

COMMON = tools_common.c tools_cell_txt.c
COMMON_O = $(COMMON:%.c=%.o)
HEADERS = tools_common.h tools_cell_txt.h

all install: $(PROG).h $(PROG).a

$(PROG).a: $(COMMON_O)
	$(AR) $(PROG).a $(COMMON_O)

$(PROG).h: $(HEADERS)
	cat $(HEADERS) > $(PROG).h

$(COMMON_O): $(COMMON) $(HEADERS)
	$(CC) $(CFLAGS) $(LIBCFLAGS) -c  $(COMMON)

clean:
	rm -f *.o $(PROG).h $(PROG).a
