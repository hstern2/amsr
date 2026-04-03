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

TARGETS := amsr/cost_grad.$(EXT) amsr/geom_embed.$(EXT)

.PHONY: all clean

all: $(TARGETS)

amsr/cost_grad.$(EXT): amsr/cost_grad.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

amsr/geom_embed.$(EXT): amsr/geom_embed.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm -f amsr/cost_grad.dylib amsr/cost_grad.so
	rm -f amsr/geom_embed.dylib amsr/geom_embed.so
