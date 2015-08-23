# Mockspec/anal/sysanal
# all necessary files are local includes; no library linking is
# required and so the library variables are commented out
FC=f95
FFLAGS=-O 
#FFLAGS=-O -Bstatic
#TLIB=/home/matrix/cwc/Programs/tlib/tool.a
#MGOLIBS=-L/usr/X11R6/lib -lmongo -lX11 -lreadline -lhistory -lcurses
#SLATEC=/home/matrix/cwc/Programs/slatec/libslatec.a -llapack
#       Insure that make knows what to do with C preprocessored fortran
.F.o: ; $(FC) $(FFLAGS) -c $*.F
#

default:: anal

anal: const.dek files.f zabs.f features.f spectra.f sysanal.com sysanal.f sysanal.h subs.f
	$(FC) $(FFLAGS) sysanal.f -o sysanal 
