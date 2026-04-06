# Build the required C extension for conformer generation.
# Usage: make          (builds for current platform)
#        make clean    (removes built libraries)

UNAME := $(shell uname -s)
CC ?= cc
CFLAGS := -O3 -fPIC
CSRC := amsr/src

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

$(TARGET): $(CSRC)/conf_util.c $(CSRC)/lbfgs.c $(CSRC)/lbfgs.h $(CSRC)/arithmetic_ansi.h
	$(CC) $(CFLAGS) -I$(CSRC) $(LDFLAGS) -o $@ $(CSRC)/conf_util.c $(CSRC)/lbfgs.c

clean:
	rm -f amsr/conf_util.dylib amsr/conf_util.so
