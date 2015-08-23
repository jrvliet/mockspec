# Mockspec/los Makefile
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

default:: specsynth
#

specsynth: specsynth.f specsynth.h const.dek convolve.f fft.f instrument.f voigt.f model.f noise.f files.f strings.f splinecode.f
	$(FC) $(OPTIONS) specsynth.f -o specsynth 
#
