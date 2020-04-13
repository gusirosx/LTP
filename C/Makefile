CC = gcc
NAME = LTP
#DEBUG
#FFLAGS=-c -fbounds-check

#RUN
CFLAGS   = -I.
DEPS     = global.h
CSOURCES = LTP_schemes.o LTP_functions.o LTP_routines.o LTP_main.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(NAME): $(CSOURCES)
	$(CC) -o $@ $^ $(CFLAGS)

run: $(NAME)
	./$(NAME)
	
clean:
	rm -f *.o
	rm -f *.exe
