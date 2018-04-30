
PROG = libezgdal

COMMON =  ezgdal.c
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



#default: $(PROG).so


#all install: $(PROG).h

#$(PROG).so: $(COMMON_O)
#	$(CC) -shared -Wl,-soname,$(PROG).so -o $(PROG).so $(COMMON_O) $(CFLAGS) 
	
#$(COMMON_O): $(COMMON) $(HEADERS)
#	$(CC) $(CFLAGS) -fPIC -c  $(COMMON)
	
#install: 
#	mv -f $(PROG).* ../../../lib
#	cp -f *.h ../../../include

#clean:
#	rm -f *.o $(PROG).a $(PROG).so
