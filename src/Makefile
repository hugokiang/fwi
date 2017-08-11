CC = gcc
OPTC = -O3
BIN = ../bin
all : propaga2d

.SUFFIXES: .c .o

.c.o:
	$(CC) $(OPTC) -c $*.c

propaga2d : propaga2d.o utils.o libfwi.o
	$(CC) $(OPTC) -o $(BIN)/propaga2d propaga2d.o utils.o libfwi.o -lm

