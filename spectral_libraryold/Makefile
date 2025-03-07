####
# Spectral test Makefile
# Generic template Makefile
#
# Builds SpectralTest statically linked for now - this is likely how Pete's dev
# process works.
#
# 'shared' target builds a shared library for testing integration with python
# GUI frontend
#
# Installs to:
# ~/meteor/bin (binaries)
# ~/meteor/include (headers)
# ~/meteor/lib
####

# Use gcc as our compiler
CC=gcc

# Default flags
CFLAGS=-fPIC -fPIC -Wall -Wno-unknown-pragmas -std=gnu99 -pthread
LDFLAGS=-pthread
LDLIBS=-lm

# Since we're still testing, set flags for debugging
CFLAGS+=-Og -ggdb


# Settings for building SpectralTest binary
SPECTRAL_BINARY=SpectralTest.so
SPECTRAL_BUILD=SpectralTestDriver.o Spectral_Functions.o System_FileFunctions.o
SPECTRAL_LDFLAGS=$(LDFLAGS)
SPECTRAL_LDLIBS=$(LDLIBS)

# Settings for building SpectralTest shared library
SHARED_LIB=SpectralTest.so.0
SHARED_SONAME=SpectralTest.so
SHARED_BUILD=Spectral_Functions.o System_FileFunctions.o
SHARED_HEADERS=Spectral_Functions.h System_FileFunctions.h StarSpectra.h
SHARED_LDFLAGS=-shared -Wl,-soname,$(SHARED_SONAME) -Wl,--no-undefined
SHARED_LDFLAGS+=$(LDFLAGS)
SHARED_LDLIBS=$(LDLIBS)

# Files for install targets
INSTALL_BIN=$(SPECTRAL_BINARY)
INSTALL_LIB=$(SHARED_LIB) $(SHARED_SONAME)
INSTALL_INC=$(SHARED_HEADERS)


#  Make rules
# ------------

default: 	spectraltest shared

spectraltest:	$(SPECTRAL_BUILD)
		$(CC) $(SPECTRAL_LDFLAGS) -o $(SPECTRAL_BINARY) $(SPECTRAL_BUILD) $(SPECTRAL_LDLIBS)

shared:		$(SHARED_BUILD)
		$(CC) $(SHARED_LDFLAGS) -o $(SHARED_LIB) $(SHARED_BUILD) $(SHARED_LDLIBS)
		ln -sf $(SHARED_LIB) $(SHARED_SONAME)

%o:		%.c
		$(CC) $(CFLAGS) -o $@ -c $<

install:	spectraltest shared
		install -m 0755 -d $(HOME)/meteor/bin
		install -m 0755 -d $(HOME)/meteor/include/spectral
		install -m 0755 -d $(HOME)/meteor/lib/spectral
		install -m 0644 -t $(HOME)/meteor/include/spectral $(INSTALL_INC)
		install -m 0755 -t $(HOME)/meteor/bin $(INSTALL_BIN)
		install -m 0755 -t $(HOME)/meteor/lib $(INSTALL_LIB)

changelog:
		git log > Changelog

tarball:	changelog
		$(eval SHORTHASH := $(shell git rev-parse --short=8 HEAD))
		mkdir -p tarballs
		git archive --prefix spectral/ --add-file Changelog -o tarballs/spectral-$(SHORTHASH).tar HEAD
		gzip tarballs/spectral-$(SHORTHASH).tar

clean:
		rm -f $(SPECTRAL_BUILD) $(SPECTRAL_BINARY)
		rm -f $(SHARED_BUILD) $(SHARED_LIB) $(SHARED_SONAME)

distclean:	clean
		git clean -fxd .
