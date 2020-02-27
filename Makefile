#
## Makefile
#
SOURCES.c= json_poly_intersect.c
$(CC)=gcc 
CFLAGS=-Wall -g -I/usr/local/include -I/usr/include -I/usr/local/jansson/include
LDLIBS= -lm -ljansson -lxml2 -L/usr/local/jansson/lib

PROGRAM=json_poly_intersect

OBJECTS= $(SOURCES.c:.c=.o)

.KEEP_STATE:

debug := CFLAGS= -g -s

all debug: $(PROGRAM)

$(PROGRAM): $(INCLUDES) $(OBJECTS)
		$(LINK.c) -o $@ $(OBJECTS) $(LDLIBS)

clean:
		rm -f $(PROGRAM) $(OBJECTS)

