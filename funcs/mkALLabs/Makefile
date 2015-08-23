# Mockspec/cullabs Makefile
# all necessary files are local includes; no library linking is
# required and so the library variables are commented out
FFLAGS=-O 
FC=f95
#FFLAGS=-O -Bstatic
#TLIB=/home/matrix/cwc/Programs/tlib/tool.a
#MGOLIBS=-L/usr/X11R6/lib -lmongo -lX11 -lreadline -lhistory -lcurses
#SLATEC=/home/matrix/cwc/Programs/slatec/libslatec.a -llapack
#       Insure that make knows what to do with C preprocessored fortran
.F.o: ; $(FC) $(FFLAGS) -c $*.F
#

default:: cullabs

cullabs: cullabs.h cullabs.f access.f files.f subs.f parze.f index.f 
	$(FC) $(FFLAGS) cullabs.f -o cullabs 
