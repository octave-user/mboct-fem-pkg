##  Copyright (C) 2018(-2020) Reinhard <octave-user@a1.net>

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

if USE_MESH_SIZE
  bin_PROGRAMS:=fem_pre_mesh_size
  fem_pre_mesh_size_SOURCES:=fem_pre_mesh_size.cc
endif

oct_OBJ:=__mboct_fem_pkg__.@OBJEXT@

oct_FILE:=__mboct_fem_pkg__.oct

bin_FILES:=$(oct_FILE)

MKOCTCXXFLAGS:=$(CXXFLAGS)
MKOCTLDFLAGS:=$(LDFLAGS)
MKOCTLIBS:=$(LIBS)

if USE_NLOPT
  MKOCTCXXFLAGS+=$(NLOPT_INC)
  MKOCTLIBS+=$(NLOPT_LIBS)
endif

if USE_OCTLIBS
  MKOCTLDFLAGS+=$(OCTLDFLAGS)
  MKOCTLIBS+=$(OCTLIBS)
endif

all-local: $(bin_FILES)

clean-local:
	rm -f $(bin_FILES) $(oct_OBJ)

distclean-local:
	rm -f oct-alt-includes.h

__mboct_fem_pkg__.@OBJEXT@: __mboct_fem_pkg__.cc
	$(MKOCTFILE) -I. $(MKOCTCXXFLAGS) -c $< -o $@

$(oct_FILE): $(oct_OBJ)
	$(MKOCTFILE) $^ $(MKOCTLDFLAGS) $(MKOCTLIBS) -o $@

.PHONY: clean-local distclean-local
