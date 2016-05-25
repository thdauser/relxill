# -*- mode: Make -*-

CFLAGS = -g -ansi -std=c99 -O2 -Wall -Wstrict-prototypes -pedantic -pg
LIBS = -L${HEADAS}/lib
LDFLAGS = -g -W -Wall $(LIBS) -lm -lcfitsio_3.37 -pg

COMPILE.c = gcc

INCLUDES = -I${HEADAS}/include

objects = test_sta.o relbase.o relmodels.o relutility.o reltable.o
headers = relbase.h  relmodels.h relutility.h reltable.h

LINK_TARGET = test_sta

.PHONY:all
all:
	make test_sta

$(LINK_TARGET): $(objects)
	gcc -o $@ $^ $(LDFLAGS) 

%.o: %.c %.h
	gcc $(INCLUDES) $(CFLAGS) -c $<

%.o: %.c
	gcc $(INCLUDES) $(CFLAGS) -c $<

clean:
	rm -f $(objects) $(LINK_TARGET) *~ gmon.out


.PHONY: valgrind, gdb
valgrind: 
	make test_sta
	valgrind --leak-check=full ./test_sta	

gdb:
	make test_sta
	gdb --args ./test_sta
