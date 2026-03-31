##############################################################################
#   PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
#   Copyright (C) 2007-2016 University of Tennessee / UT-Battelle
#
#   Code Author:  Doug Hyatt
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

SHELL   = /bin/sh
CC      = gcc

CFLAGS  += -pedantic -Wall -O3 -DSUPPORT_GZIP_COMPRESSED
LFLAGS = -lm $(LDFLAGS) -lz

TARGET  = prodigal
INSTALLDIR  = /usr/local/bin

# Source file groups
CORE_SOURCES = bitmap.c dprog.c gene.c metagenomic.c node.c sequence.c training.c
API_SOURCES  = prodigal_api.c
CLI_SOURCE   = main.c
TEST_SOURCE  = test_api.c

HEADERS = $(shell echo *.h)

CORE_OBJS = $(CORE_SOURCES:.c=.o)
API_OBJ   = $(API_SOURCES:.c=.o)
LIB_OBJS  = $(CORE_OBJS) $(API_OBJ)
CLI_OBJ   = $(CLI_SOURCE:.c=.o)

# Default: build CLI binary
all: $(TARGET)

# Core and API objects
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

# Static library (no main, no zlib dependency in library itself)
libprodigal.a: $(LIB_OBJS)
	ar rcs $@ $^

# CLI binary: link main.o with static library
$(TARGET): $(CLI_OBJ) libprodigal.a
	$(CC) $(CFLAGS) -o $@ $(CLI_OBJ) -L. -lprodigal $(LFLAGS)

# PIC objects for shared library
%.pic.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

LIB_PIC_OBJS = $(CORE_SOURCES:.c=.pic.o) $(API_SOURCES:.c=.pic.o)

# Shared library
libprodigal.so: $(LIB_PIC_OBJS)
	$(CC) -shared -o $@ $^ -lm

# Test runner
test_api: $(TEST_SOURCE) libprodigal.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lprodigal $(LFLAGS)

test: test_api
	./test_api

install: $(TARGET)
	install -d -m 0755 $(INSTALLDIR)
	install -m 0755 $(TARGET) $(INSTALLDIR)

uninstall:
	-rm $(INSTALLDIR)/$(TARGET)

clean:
	-rm -f *.o *.pic.o

distclean: clean
	-rm -f $(TARGET) libprodigal.a libprodigal.so test_api

.PHONY: all install uninstall clean distclean test
