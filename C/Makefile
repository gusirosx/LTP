PROGRAM = LTP
CC = gcc
#DEBUG
#FFLAGS=-c -fbounds-check

#RUN
CFLAGS   = -I.
#CFLAGS=-I. -O3
DEPS     = global.h
CSOURCES = LTP_schemes.o LTP_functions.o LTP_routines.o LTP_main.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(PROGRAM): $(CSOURCES)
	$(CC) -o $@ $^ $(CFLAGS)

run: $(PROGRAM)
	./$(PROGRAM)
	
clean:
	rm -f *.o
	rm -f *.exe
