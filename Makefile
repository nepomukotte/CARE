# Makefile for CARE
#
#
# 
#



ARCH := $(shell uname)

# linux flags
ifeq ($(ARCH),Linux)
CXX           = g++ 
CXXFLAGS      =  -g -O3 -Wall -fPIC -fno-strict-aliasing  -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -D_LARGEFILE64_SOURCE
LD            = g++
LDFLAGS       = -g
# LDFLAGS       =  -pg -O
SOFLAGS       = -shared
endif

# Apple OS X flags
ifeq ($(ARCH),Darwin)
CXX           = g++ 
CXXFLAGS      = -g -O3 -Wall -fPIC  -fno-strict-aliasing
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

OutPutOpt     = -o


# get vbf details
CPPFLAGS   := $(shell vbfConfig --cppflags)
CXXFLAGS   := $(shell vbfConfig --cxxflags)
CLLFLAGS   := $(shell vbfConfig --ldflags)
VBFLIBS    := $(shell vbfConfig --libs)

# Root

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLDFLAGS  := $(shell root-config --ldflags)

CLLFLAGS += $(ROOTLDFLAGS)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTGLIBS)
LIBS         += -lMinuit
LIBS         += $(VBFLIBS)

INCLUDEFLAGS  = -I. -I./inc/
CXXFLAGS     += $(INCLUDEFLAGS)

ALLFLAGS = $(CXXFLAGS) $(CPPFLAGS) -Wall


# rule for any compiling any .cpp file
%.o : %.cpp  
	@printf "Compiling $< ... "
	@g++ $(ALLFLAGS) -c $<
	@echo "Done"

CameraAndReadout: GOrderedGrid.o  GOrderedGridSearch.o VG_writeVBF.o VATime.o CameraAndReadout.o TelescopeData.o TraceGenerator.o TriggerTelescopeBase.o TriggerTelescopeVERITAS.o TriggerTelescopeSPB2.o TriggerTelescopeTriDem.o ArrayTrigger.o ReadConfig.o FADC.o Display.o
	        $(LD)  $(CLLFLAGS) $^ $(LIBS)  $(LDFLAGS) $(OutPutOpt) $@ -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS
	        @echo "$@ done"

VATime.o: VATime.cpp
	g++ $(ALLFLAGS) VATime.cpp -c -o VATime.o -DNOROOT -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS

clean:	
	rm -f *.o

.SUFFIXES: .o


depend:
	$(DEPEND)

Makefile.depend:
	$(DEPEND)

DEPEND=echo > Makefile.depend0 &&\
	makedepend -s "\#DEPEND LIST DONT DELETE" -- $(INCLUDEFLAGS) \
		-Y --  *.cpp *.c \
		-f Makefile.depend0 > /dev/null 2>&1 &&\
	sed "s/^[a-zA-Z0-9]*\///" Makefile.depend0 > Makefile.depend &&\
	rm -f Makefile.depend0

include Makefile.depend
