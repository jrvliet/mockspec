this = .
PROG = ana.R.exe

SRCS =  analysis_Rec.F ART_Control.F analysis_ART.F \
	 analysis_Auxiliary.F analysis_SubR.F \
	ART_Init.F

OBJS =  analysis_Rec.o ART_Control.o analysis_ART.o \
	 analysis_Auxiliary.o analysis_SubR.o \
	ART_Init.o

includeall = a_tree.h a_control.h a_def.h analysis_Rec.h a_setup.h

# Intel compilers:
FC = ifort
#FC = /home/local/intel/intel_fc_80/bin/ifort
#FFLAGS =  -O3 -axW -ftz -unroll -openmp -mcmodel=medium -i-dynamic -xW -convert big_endian 
#LDFLAGS = -O3 -axW -ftz -unroll -openmp -mcmodel=medium -i-dynamic -xW -convert big_endian

# What was the right ones
FFLAGS = -O3      -convert big_endian -mcmodel=medium  -i-dynamic -xW -axW
LDFLAGS= -O3      -convert big_endian -mcmodel=medium  -i-dynamic -xW -axW


all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $(this)/$(PROG) $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS)

analysis_Rec.o: analysis_Rec.F $(includeall)
	$(FC) $(FFLAGS) -c analysis_Rec.F
analysis_SubR.o: analysis_SubR.F $(includeall)
	$(FC) $(FFLAGS) -c analysis_SubR.F
ART_Control.o: ART_Control.F a_tree.h a_control.h a_numbers.h a_constant.h $(includeall)
	$(FC) $(FFLAGS) -c ART_Control.F $(LIB)
analysis_ART.o: analysis_ART.F $(includeall)
	$(FC) $(FFLAGS) -c analysis_ART.F
analysis_Auxiliary.o: analysis_Auxiliary.F $(includeall)
	$(FC) $(FFLAGS) -c analysis_Auxiliary.F
ART_Init.o: ART_Init.F a_tree.h a_control.h a_numbers.h a_constant.h $(includeall)
	$(FC) $(FFLAGS) -c ART_Init.F $(LIB)
