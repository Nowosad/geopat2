
PROG = libsml

COMMON =  smllib.c smltxt.c
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
	mv -f $(PROG).* ../../../lib
	cp -f *.h ../../../include

clean:
	rm -f *.o $(PROG).h $(PROG).a

#dynamic: smllib.c smltxt.c sml.h
#	$(CC) -c -fPIC smllib.c -o smllib.o
#	$(CC) -c -fPIC smltxt.c -o smltxt.o
#	$(CC) -shared -Wl,-soname,$(PROG).so -o $(PROG).so smllib.o smltxt.o

#static: sml.h smllib.c
#	$(CC) -fdata-sections -ffunction-sections -Wl,--gc-sections -c smllib.c -o smllib.o
#	$(AR) $(PROG).a smllib.o
	
#$(PROG).so: $(COMMON_O)
#	$(CC) -shared -Wl,-soname,$(PROG).so -o $(PROG).so $(COMMON_O) $(CFLAGS) 
	
#$(COMMON_O): $(COMMON) $(HEADERS)
#	$(CC) $(CFLAGS) -fPIC -c  $(COMMON)
	
#all install: $(PROG).a

#all install: 
#	cp -f $(PROG).* $(PREFIX)/lib
#	cp -f sml.h $(PREFIX)/include/sml.h

#clean:
#	rm -f *.o $(PROG).a $(PROG).so
