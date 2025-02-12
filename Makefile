TARGET := svtrek
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

# directories
CURRENT_DIR := $(shell pwd)
BIN_DIR := $(CURRENT_DIR)/bin

# compiler
GXX := g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++11
TIME := /usr/bin/time -v

# object files that need lcptools
HTSLIB_CXXFLAGS := -I$(CURRENT_DIR)/htslib/include
HTSLIB_LDFLAGS := -L$(CURRENT_DIR)/htslib/lib -lhts -Wl,-rpath,$(CURRENT_DIR)/htslib/lib -pthread

$(TARGET): $(OBJS)
	$(GXX) $(CXXFLAGS) -o $@ $^ $(HTSLIB_LDFLAGS)
	rm *.o

%.o: %.cpp
	$(GXX) $(CXXFLAGS) $(HTSLIB_CXXFLAGS) -c $< -o $@

clean:
	rm -rf *.o svtrek
