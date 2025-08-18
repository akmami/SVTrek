TARGET := svtrek
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)

# directories
CURRENT_DIR := $(shell pwd)
BIN_DIR := $(CURRENT_DIR)/bin

# compiler
GXX := gcc
CXXFLAGS = -Wall -Wextra -O3
TIME := /usr/bin/time -v

# object files that need lcptools
HTSLIB_CXXFLAGS := -I$(CURRENT_DIR)/htslib/include
HTSLIB_LDFLAGS := -L$(CURRENT_DIR)/htslib/lib -lhts -Wl,-rpath,$(CURRENT_DIR)/htslib/lib -pthread -lz

$(TARGET): $(OBJS)
	$(GXX) $(CXXFLAGS) -o $@ $^ $(HTSLIB_LDFLAGS)
	rm *.o

%.o: %.c
	$(GXX) $(CXXFLAGS) $(HTSLIB_CXXFLAGS) -c $< -o $@

install:
	@echo "Installing htslib"
	cd htslib && \
	autoreconf -i && \
	./configure && \
	make && \
	make prefix=$(CURRENT_DIR)/htslib install && \
	cd ../abPOA && \
	make


clean:
	rm -rf *.o svtrek
