# Build the required C extension for conformer generation.
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

TARGET := amsr/conf_util.$(EXT)

.PHONY: all clean

all: $(TARGET)

$(TARGET): amsr/conf_util.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm -f amsr/conf_util.dylib amsr/conf_util.so
