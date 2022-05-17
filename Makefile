CC=gcc
CFLAGS =  -O3 -funroll-loops -g -I htslib
LDFLAGS = htslib/libhts.a -lz -lm -lpthread -llzma -lbz2
SOURCES = sveldt.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = sveldt

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@
