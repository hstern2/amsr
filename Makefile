# Build the required C extension for cost/gradient evaluation.
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

TARGET := amsr/_cost_grad_c.$(EXT)

.PHONY: all clean

all: $(TARGET)

$(TARGET): amsr/_cost_grad_c.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm -f amsr/_cost_grad_c.dylib amsr/_cost_grad_c.so
