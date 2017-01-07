# -*- mode: Make -*-

CFLAGS = -g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic -O2
LDFLAGS = -g -W -Wall $(LIBS) -lm -lcfitsio

LIBS = -L${HEADAS}/lib

COMPILE.c = gcc

INCLUDES = -I${HEADAS}/include

objects = test_sta.o relbase.o relmodels.o relutility.o reltable.o rellp.o
headers = relbase.h  relmodels.h relutility.h reltable.h rellp.h common.h

LINK_TARGET = test_sta

.PHONY:all
all:
	make test_sta

$(LINK_TARGET): $(objects)
	echo $(CFLAGS)
	gcc -o $@ $^ $(LDFLAGS) 

%.o: %.c %.h
	gcc $(INCLUDES) $(CFLAGS) -c $<

%.o: %.c
	gcc $(INCLUDES) $(CFLAGS) -c $<

clean:
	rm -f $(objects) $(LINK_TARGET) *~ gmon.out


.PHONY: valgrind, gdb
valgrind:
	make clean
	make CFLAGS="-g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic" test_sta
	valgrind --tool=memcheck --leak-check=full ./test_sta	

gdb:
	make clean
	make CFLAGS="-g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic" test_sta
	gdb --args ./test_sta

ddd:
	echo "exit" | make gdb 
	ddd test_sta &

gprof:
	make clean
	make CFLAGS="$(CFLAGS) -pg" LDFLAGS="$(LDFLAGS) -pg" test_sta 
	./test_sta
	gprof -p test_sta
