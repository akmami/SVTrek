CC=g++
CFLAGS =  -O3 -funroll-loops -g -I htslib
LDFLAGS = htslib/libhts.a -lz -lm -lpthread -llzma -lbz2 -lcurl
SOURCES = sveldt.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = sveldt

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@
