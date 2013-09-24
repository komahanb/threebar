# Copyright (C) 2003, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# $Id: Makefile.in 1875 2010-12-28 23:32:54Z andreasw $

##########################################################################
#    You can modify this erpxample makefile to fit for your own program.   #
#    Usually, you only need to change the five CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable

EXE = threebar	

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS = threebar.o

# CHANGEME: Additional libraries
ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

# Fortran Compiler options
F77 = ifort

# Fotran Compiler options
FFLAGS =   -r8 -check -debug all -traceback -ftrapuv -check uninit
# additional Fortran Compiler options for linking
F77LINKFLAGS =  -Wl,--rpath -Wl, /home/komahan/Dropbox/Thesis/Program/Markus/Ipopt-3.10.0/lib


# Linker flags
LIBS = `PKG_CONFIG_PATH=/home/komahan/Dropbox/Thesis/Program/Markus/Ipopt-3.10.0/lib/pkgconfig:/home/komahan/Dropbox/Thesis/Program/Markus/Ipopt-3.10.0/share/pkgconfig: /usr/bin/pkg-config --libs ipopt` -lstdc++ -lm


all: $(EXE)

.SUFFIXES: .f90 .o

$(EXE): $(OBJS)
	$(F77) $(F77LINKFLAGS) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(EXE) $(OBJS) IPOPT.OUT *~

.f90.o:
	$(F77) $(FFLAGS) -c -o $@ $<

