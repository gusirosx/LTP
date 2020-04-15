PROGRAM = LTP
FC   = gfortran

#DEBUG
#FFLAGS=-c -fbounds-check

#RUN
FFLAGS=-c
#FFLAGS=-c -O3
#FFLAGS=-c -O4
#FFLAGS=-c -O5

FSOURCES=LTP_main.o LTP_global.o LTP_functions.o LTP_schemes.o LTP_routines.o

$(PROGRAM): $(FSOURCES)
	$(FC)  $(FSOURCES) -o $(PROGRAM)

LTP_global.mod: LTP_global.f90 LTP_global.o 
	$(FC) $(FFLAGS) LTP_global.f90

LTP_global.o: LTP_global.f90
	$(FC) $(FFLAGS) LTP_global.f90

LTP_schemes.mod: LTP_schemes.f90 LTP_schemes.o LTP_global.mod 
	$(FC) $(FFLAGS) LTP_schemes.f90

LTP_schemes.o: LTP_schemes.f90 LTP_global.mod
	$(FC) $(FFLAGS) LTP_schemes.f90

LTP_functions.mod: LTP_functions.f90 LTP_functions.o LTP_global.mod LTP_schemes.mod
	$(FC) $(FFLAGS) LTP_functions.f90

LTP_functions.o: LTP_functions.f90 LTP_global.mod LTP_schemes.mod
	$(FC) $(FFLAGS) LTP_functions.f90
	
LTP_routines.mod: LTP_routines.f90 LTP_routines.o LTP_global.mod LTP_schemes.mod LTP_functions.mod
	$(FC) $(FFLAGS) LTP_routines.f90

LTP_routines.o: LTP_routines.f90 LTP_global.mod LTP_schemes.mod LTP_functions.mod
	$(FC) $(FFLAGS) LTP_routines.f90

LTP_main.o: LTP_main.f90  LTP_routines.mod
	$(FC) $(FFLAGS) LTP_main.f90

run: $(PROGRAM)
	./$(PROGRAM)

clean:
	rm -f *.f90~
	rm -f *.mod
	rm -f *.o
	rm -f *.exe
