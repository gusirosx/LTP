PROGRAM = LTP
FC   = gfortran

#DEBUG
#FFLAGS=-c -fbounds-check

#RUN
FFLAGS=-c
#FFLAGS=-c -g -Wall -Werror -Wconversion -Wno-unused-dummy-argument -fbacktrace -fcheck=all
#FFLAGS=-c -O3
#FFLAGS=-c -O4
#FFLAGS=-c -O5

#FOLDERS
OUT = ./output
BIN = ./bin
SRC = ./sources

#SOURCES
MAIN = $(SRC)/LTP_main.f90
FSOURCES = $(BIN)/LTP_global.o $(BIN)/LTP_schemes.o $(BIN)/LTP_functions.o $(BIN)/LTP_routines.o

all: create $(PROGRAM)

$(PROGRAM): $(MAIN) $(FSOURCES)
	$(FC)  $(FSOURCES) $(MAIN) -I $(BIN) -o $(BIN)/$(PROGRAM)
	@echo ' '
	@echo '--- Finished Building Program: $@ ---'

$(BIN)/%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) $< -J $(BIN) -o $@

run: create $(PROGRAM)
	./bin/$(PROGRAM)
	
clean:
	rm -f $(BIN)/*
	rm -f $(OUT)/*
	
delete:
	rm -r $(BIN)
	rm -r $(OUT)
	
create:
	@mkdir -p $(OUT) $(BIN)

cleanout:
	rm -f $(OUT)/*
	
#Backup
date=` data +%F `
version=` date +%d.%m.%y_%H.%M`

output      = ../$(PROGRAM)_$(version).tar.gz
saveindir   = ../$(PROGRAM)_INPUT_$(version).tar.gz
casedir     = ../$(PROGRAM)_CASE_$(version).tar.gz
files       = ./sources/*.f90 ./input/*.ltp ./makefile
files_input = ./input/*.ltp
files_case  = ./input/*.ltp ./output/*
backup:
	@tar -czvf $(output) $(files)
	@echo '--- All program files have been saved ---'
	
saveinput: 
	@tar -czvf $(saveindir) $(files_input)
	@echo '--- All input files have been saved ---'
	
savecase: 
	@tar -czvf $(casedir) $(files_case)
	@echo '--- All case files have been saved ---'
	
help: 
	@echo '_________________________________________________________________________'
	@echo '                                                                         '
	@echo '                           Make Assistant Tool                           '
	@echo '_________________________________________________________________________'
	@echo '                               Make Flags                                '
	@echo '         make  < flag >                                                  '
	@echo '               run           : build and execute'
	@echo '               clean         : clears compilation files'
	@echo '               delete        : delete directories created'
	@echo '               cleanout      : clean results directory'
	@echo '               backup        : performs a code backup'
	@echo '               saveinput     : backup all input files'
	@echo '               savecase      : backup all case files'
	@echo '_________________________________________________________________________'

.PHONY: clean backup delete run
