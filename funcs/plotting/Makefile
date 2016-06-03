FFLAGS=-O 
#FFLAGS=-O -Bstatic
TLIB=./tool.a
OPTIONS = -O
FPP=f
FC=f95
.SUFFIXES:
.SUFFIXES: .out .o .c .e .r .F .f .y .l .s .p
#       Insure that make knows what to do with C preprocessored fortran
.F.o: ; $(FC) $(FFLAGS) -c $*.F
#


bindata-logfreq: bindata-logfreq.f
	$(FC) $(OPTIONS) bindata-logfreq.f -o bindata-logfreq $(TLIB)

