# -*- mode: Make -*-

CFLAGS = -g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic -O3
LDFLAGS = -g -W -Wall $(LIBS) -lm -lcfitsio


LIBS = -L${HEADAS}/lib

COMPILE.c = gcc

INCLUDES = -I${HEADAS}/include

objects = test_sta.o relbase.o relmodels.o relutility.o reltable.o rellp.o xilltable.o 
headers = relbase.h  relmodels.h relutility.h reltable.h rellp.h common.h xilltable.h
sourcefiles = relbase.c  relmodels.c relutility.c reltable.c rellp.c xilltable.c

model_dir = ./build/
model_files = $(headers) $(sourcefiles) modelfiles/lmodel_relxill.dat modelfiles/compile_relxill.csh modelfiles/README.txt

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
	rm -f $(objects) $(LINK_TARGET) *~ gmon.out test*.dat
	rm -rf $(model_dir)

MODEL_VERSION = x.y.z
MODEL_TAR_NAME = relxill_model_v$(MODEL_VERSION).tgz


.PHONY: model
model:
	mkdir -p $(model_dir)
	rm -f $(model_dir)/*
	cp -v $(model_files) $(model_dir)

	make $(LINK_TARGET)
	$(eval MODEL_VERSION := $(shell ./test_sta version))
	$(eval MODEL_TAR_NAME := relxill_model_v$(MODEL_VERSION).tgz)

	cd $(model_dir) && tar cfvz $(MODEL_TAR_NAME) *
	cd $(model_dir) && ./compile_relxill.csh && echo 'require("xspec"); load_xspec_local_models("."); fit_fun("relxill"); () = eval_fun(1,2); exit; ' | isis -v 
	cp $(model_dir)/$(MODEL_TAR_NAME) .
	rm -f $(model_dir)/*.c $(model_dir)/*.h
	@echo "\n  --> Built model  *** $(MODEL_TAR_NAME) *** \n"



.PHONY: valgrind, gdb
valgrind:
	make clean
	make CFLAGS="-g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic" test_sta
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./test_sta relxill

valgrind-relxilllp:
	make clean
	make CFLAGS="-g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic" test_sta
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./test_sta relxilllp

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
