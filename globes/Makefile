# GLoBES -- General LOng Baseline Experiment Simulator
# (C) 2002 - 2007,  The GLoBES Team
#
# GLoBES is mainly intended for academic purposes. Proper
# credit must be given if you use GLoBES or parts of it. Please
# read the section 'Credit' in the README file.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# The configure script should have replaced anything as required to
# obtain a working Makefile, which may serve as a template for own
# applications.
#
# This simple Makefile is for the GLoBES examples
#
# Compile example files with ``make eventrate'' etc.
#
# This Makefile assumes dynamic libraries, installed to either
# the default prefix /usr/local/ or to a user-defined directory 
# called ${prefix}.
#
# For linking against a specific version of GLoBES, libglobes.so can be 
# replaced by the respective library, such as libglobes.so.0.0.1

#UNAME := $(shell uname)
#ifneq ($(UNAME),Darwin)
#CFLAGS := $(CFLAGS) -static
#endif


prefix = /Users/cvson/meOffline/globes/rebuiltglobes3218
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

local_CFLAGS = -g -O4
local_CFLAGS += $(shell root-config --cflags)
#local_CFLAGS += -D_REENTRANT -stdlib=libc++ -I/Users/cvson/opt/anaconda3/envs/my_root5_env/include/root 
#local_CFLAGS += -D_REENTRANT -stdlib=libc++ -I/Users/cvson/meOffline_nosync/root6/root_reinstall/include
#local_CFLAGS +=  -I/Library/Developer/CommandLineTools/usr/include/c++/v1

INCFLAGS:=$(shell $(globesconf) --include)
local_LDFLAGS:=$(shell $(globesconf) --libs) 
local_LDFLAGS +=$(shell root-config --libs) 
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)



BIN = glbProb eventrateall test_cp_null sensi_th23dm31 sensi_th23dcp sensi_th23cosdcp sensifut_th23cosdcp example7
OBJ = glbProb.o eventrateall.o test_cp_null.o sensi_th23dm31.o sensi_th23dcp.o sensi_th23cosdcp.o sensifut_th23cosdcp.o example7.o

all: $(BIN)

glbProb: glbProb.o
	g++ glbProb.o -o glbProb  $(LDFLAGS) $(local_LDFLAGS)

eventrateall: eventrateall.o
	g++ eventrateall.o  -o eventrateall  $(LDFLAGS) $(local_LDFLAGS)

test_cp_null: test_cp_null.o
	g++ test_cp_null.o -o test_cp_null  $(LDFLAGS) $(local_LDFLAGS)

sensi_th23dm31: sensi_th23dm31.o
	g++ sensi_th23dm31.o -o sensi_th23dm31  $(LDFLAGS) $(local_LDFLAGS)

sensi_th23dcp: sensi_th23dcp.o
	g++ sensi_th23dcp.o -o sensi_th23dcp  $(LDFLAGS) $(local_LDFLAGS)
sensi_th23cosdcp: sensi_th23cosdcp.o
	g++ sensi_th23cosdcp.o -o sensi_th23cosdcp  $(LDFLAGS) $(local_LDFLAGS)
sensifut_th23cosdcp: sensifut_th23cosdcp.o
	g++ sensifut_th23cosdcp.o -o sensifut_th23cosdcp  $(LDFLAGS) $(local_LDFLAGS)
example7: example7.o
	g++ example7.o -o example7  $(LDFLAGS) $(local_LDFLAGS)
%.o : %.c
	g++ $(CFLAGS) $(local_CFLAGS) -c $< $(INCFLAGS)
.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)
