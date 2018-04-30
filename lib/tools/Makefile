
PROG = libtools

COMMON = tools_common.c tools_cell_txt.c
COMMON_O = $(COMMON:%.c=%.o)
HEADERS = $(shell ls *.h)


all install: $(PROG).h $(PROG).a

$(PROG).a: $(COMMON_O)
	$(AR) $(PROG).a $(COMMON_O)

$(PROG).h: $(HEADERS)
	cat $(HEADERS) > $(PROG).h

$(COMMON_O): $(COMMON) $(HEADERS)
	$(CC) $(CFLAGS) $(LIBCFLAGS) -c  $(COMMON)

install: $(PROG).a
	#mv -f $(PROG).* ../../../lib
	cp -f *.h ../../../include

clean:
	rm -f *.o $(PROG).h $(PROG).a
