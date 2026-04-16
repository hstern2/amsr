# Build the required C extension for conformer generation.
# Usage: make          (builds for current platform)
#        make clean    (removes built libraries)

UNAME := $(shell uname -s)
CC ?= cc
CFLAGS := -O3 -fPIC
LDFLAGS := -shared -lm
SRC := amsr/src

ifeq ($(UNAME),Darwin)
  EXT := dylib
else
  EXT := so
endif

TARGET := $(SRC)/conf_util.$(EXT)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)/conf_util.c $(SRC)/lbfgs.c $(SRC)/lbfgs.h $(SRC)/arithmetic_ansi.h
	$(CC) $(CFLAGS) -I$(SRC) $(LDFLAGS) -o $@ $(SRC)/conf_util.c $(SRC)/lbfgs.c

clean:
	rm -f $(TARGET)
