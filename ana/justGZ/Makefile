this = .
PROG = hydroANA-justGZ.x
O = .
SRCS = ART_Auxiliary.F ART_Control.F ART_Conversions.F \
	ART_Init.F ART_IO.F ART_Tools.F analysis_Rec.F analysis_SubR.F
OBJS = $(O)/ART_Auxiliary.o $(O)/ART_Control.o $(O)/ART_Conversions.o \
	$(O)/ART_Init.o $(O)/ART_IO.o $(O)/ART_Tools.o \
	$(O)/analysis_Rec.o $(O)/analysis_SubR.o
includeall = a_control.h a_def.h a_numbers.h a_setup.h a_tree.h

# Intel compilers:
FC = ifort
FFLAGS = -O2 -convert big_endian -mcmodel=medium -i-dynamic
LDFLAGS= -O2 -convert big_endian -mcmodel=medium -i-dynamic
## Blue Waters:
#FC = ftn
#FFLAGS  = -O2 -fopenmp -dynamic -mcmodel=medium -shared-intel -convert big_endian
#LDFLAGS = -O2 -fopenmp -dynamic -mcmodel=medium -shared-intel -convert big_endian


all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $(this)/$(PROG) $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS)

$(O)/ART_Auxiliary.o: ART_Auxiliary.F
	$(FC) $(FFLAGS) -c ART_Auxiliary.F
$(O)/ART_Control.o: ART_Control.F a_control.h a_def.h a_tree.h
	$(FC) $(FFLAGS) -c ART_Control.F
$(O)/ART_Conversions.o: ART_Conversions.F $(includeall)
	$(FC) $(FFLAGS) -c ART_Conversions.F
$(O)/ART_Init.o: ART_Init.F $(includeall)
	$(FC) $(FFLAGS) -c ART_Init.F
$(O)/ART_IO.o: ART_IO.F $(includeall)
	$(FC) $(FFLAGS) -c ART_IO.F
$(O)/ART_Tools.o: ART_Tools.F a_control.h a_setup.h a_tree.h
	$(FC) $(FFLAGS) -c ART_Tools.F
$(O)/analysis_Rec.o: analysis_Rec.F a_control.h a_def.h a_tree.h
	$(FC) $(FFLAGS) -c analysis_Rec.F
$(O)/analysis_SubR.o: analysis_SubR.F a_control.h a_def.h a_tree.h
	$(FC) $(FFLAGS) -c analysis_SubR.F

