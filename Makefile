# Build the optional C extension for faster cost/gradient evaluation.
# Usage: make          (builds for current platform)
#        make clean    (removes built libraries)

UNAME := $(shell uname -s)
CC ?= cc
CFLAGS := -O3 -fPIC

ifeq ($(UNAME),Darwin)
  EXT := dylib
  LDFLAGS := -shared -lm
else
  EXT := so
  LDFLAGS := -shared -lm
endif

TARGET := amsr/cost_grad.$(EXT)

.PHONY: all clean

all: $(TARGET)

$(TARGET): amsr/cost_grad.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm -f amsr/cost_grad.dylib amsr/cost_grad.so
