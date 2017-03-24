# -*- mode: Make -*-

CFLAGS = -g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic -O3
LDFLAGS = -g -W -Wall $(LIBS) -lm -lcfitsio

MODEL_VERSION = 0.1.4dev
MODEL_TAR_NAME = relxill_model_v$(MODEL_VERSION).tgz
LIBS = -L${HEADAS}/lib

COMPILE.c = gcc

INCLUDES = -I${HEADAS}/include

objects = test_sta.o relbase.o relmodels.o relutility.o reltable.o rellp.o xilltable.o 
headers = relbase.h  relmodels.h relutility.h reltable.h rellp.h common.h xilltable.h
sourcefiles = relbase.c  relmodels.c relutility.c reltable.c rellp.c xilltable.c

model_dir = ./build/
# add_model_files = modelfiles
model_files = $(headers) $(sourcefiles) modelfiles/lmodel_relxill.dat modelfiles/compile_relxill.csh modelfiles/README.txt

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
	rm -f $(objects) $(LINK_TARGET) *~ gmon.out test*.dat
	rm -rf $(model_dir)


.PHONY: model
model:
	mkdir -p $(model_dir)
	rm -f $(model_dir)/*
	cp -v $(model_files) $(model_dir)
	cd $(model_dir) && tar cfvz $(MODEL_TAR_NAME) *
	cd $(model_dir) && ./compile_relxill.csh && echo 'load_xspec_local_models("."); fit_fun("relxill"); () = eval_fun(1,2); exit; ' | isis
	cp $(model_dir)/$(MODEL_TAR_NAME) .
#	rm -f $(model_dir)/*.c $(model_dir)/*.h 


.PHONY: valgrind, gdb
valgrind:
	make clean
	make CFLAGS="-g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic" test_sta
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./test_sta relxilllp

valgrind-rellinelp:
	make clean
	make CFLAGS="-g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic" test_sta
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./test_sta rellinelp

gdb:
	make clean
	make CFLAGS="-g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic" test_sta
	gdb --args ./test_sta relxilllp

ddd:
	echo "exit" | make gdb 
	ddd test_sta &

gprof:
	make clean
	make CFLAGS="$(CFLAGS) -pg" LDFLAGS="$(LDFLAGS) -pg" test_sta 
	./test_sta relxilllp 100
	gprof -p test_sta
