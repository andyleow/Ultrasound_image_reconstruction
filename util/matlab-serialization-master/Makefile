# Makefile for mex files.
#
# Example:
#
#    make MATLABDIR=/Applications/MATLAB_R2013a.app
#
MATLABDIR   ?= /usr/local/matlab
MEX         ?= $(MATLABDIR)/bin/mex
MEXEXT      ?= $(shell $(MATLABDIR)/bin/mexext)
RM          ?= rm
SOURCES     := $(wildcard *.cc)
TARGETS     := $(SOURCES:.cc=.$(MEXEXT))

all: $(TARGETS)

%.$(MEXEXT): %.cc
	$(MEX) $<

clean:
	$(RM) -rf *.$(MEXEXT)
