FFLAGS=-O -C -w 
#FFLAGS=-O -Bstatic
#SLATEC=/home/matrix/cwc/Programs/slatec/libslatec.a -llapack
SLATEC=/home/matrix/cwc/Swinburne/Swin-Programs/Programs/slatec/libslatec.a #-llapack# -dlamch
LAPACK=/home/matrix/cwc/Swinburne/Swin-Programs/Programs/lapack/CLAPACK/lapack_LINUX.a /home/matrix/cwc/Swinburne/Swin-Programs/Programs/lapack/CLAPACK/F2CLIBS/libF77.a  /home/matrix/cwc/Swinburne/Swin-Programs/Programs/lapack/CLAPACK/F2CLIBS/libI77.a
#       Insure that make knows what to do with C preprocessored fortran
.F.o: ; f95 $(FFLAGS) -c $*.F
#

default:: rates

rates: files-auger.f files-coll.f files-photo.f files-recomb.f funcs-auger.f funcs-coll.f funcs-photo.f funcs-recomb.f getcube.f kingdon-mod.f phfit-mod.f rates-alphas.f rates-phxsecs.f rates-comm.f rates.f rates-files.f rates-sed.f rates-solve4.f rates-thick.f rates-subs.f rates-times.f com-auger.com com-coll.com com-modes.com com-photo.com com-recomb.com getcube.com rates.com getcube.h rates.h

	f95 $(FFLAGS) rates.f -o rates2 $(SLATEC) $(LAPACK)
