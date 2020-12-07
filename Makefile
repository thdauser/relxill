# -*- mode: Make -*-

RELXILL_SOURCE_DIR = ./

CFLAGS = -g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic -O3 -DRELXILL_SOURCE_DIR='"${RELXILL_SOURCE_DIR}"'

LDFLAGS = -g -W -Wall $(LIBS) -lm -lcfitsio


LIBS = -L${HEADAS}/lib

COMPILE.c = gcc

INCLUDES = -I/usr/include -I${HEADAS}/include 

objects = test_sta.o relbase.o relmodels.o relutility.o reltable.o rellp.o xilltable.o donthcomp.o relcache.o test_relxill.o relphysics.o test_rellp.o test_std_functions.o
headers = relbase.h  relmodels.h relutility.h reltable.h rellp.h common.h test_relxill.h xilltable.h relcache.h relphysics.h test_rellp.h test_std_functions.h
sourcefiles = relbase.c  relmodels.c relutility.c reltable.c rellp.c test_relxill.c xilltable.c donthcomp.c relcache.c test_xilltab.c relphysics.c test_rellp.c test_std_functions.c

model_dir = ./build/
model_files = $(headers) $(sourcefiles) modelfiles/lmodel_relxill.dat modelfiles/compile_relxill.sh modelfiles/README.txt modelfiles/CHANGELOG.txt

unpublished_model_files = relreturn.c relreturn.h relreturn_corona.c  relreturn_corona.h  relreturn_datastruct.c  relreturn_datastruct.h  relreturn_table.c  relreturn_table.h

LINK_TARGET = test_sta

.PHONY:all
all:
	echo $(objects)
	make test_sta


.PHONY: test
test:
	cd test/e2e && make test


$(LINK_TARGET): $(objects)
	gcc -o $@ $^ $(LDFLAGS) 

%.o: %.c %.h
	gcc $(INCLUDES) $(CFLAGS) -c $<

%.o: %.c
	gcc $(INCLUDES) $(CFLAGS) -c $<

clean:
	rm -f $(objects) $(LINK_TARGET) *~ gmon.out test*.dat *.log
	rm -rf $(model_dir)
	rm -f relxill_model_v*.tgz
	rm -f debug-*fits testrr-*.fits

MODEL_VERSION = undef
#MODEL_TAR_NAME = relxill_model_v$(MODEL_VERSION).tgz


.PHONY: compilemodel
compilemodel: test_sta

	$(eval MODEL_TAR_NAME := relxill_model_v$(MODEL_VERSION)$(DEV).tgz)

	cd $(model_dir) && tar cfvz $(MODEL_TAR_NAME) *

	cd $(model_dir) && ./compile_relxill.sh && echo 'require("xspec"); load_xspec_local_models("./librelxill.so"); fit_fun("relxill"); () = eval_fun(1,2); exit; ' | isis -v
	cp $(model_dir)/$(MODEL_TAR_NAME) .
	rm -f $(model_dir)/*.c $(model_dir)/*.h
	@echo "\n  --> Built model  *** $(MODEL_TAR_NAME) *** \n"


.PHONY: model
model: test_sta
	mkdir -p $(model_dir)
	rm -f $(model_dir)/*
	cp -v $(model_files) $(model_dir)

	$(eval MODEL_VERSION := $(shell ./test_sta version))
	make compilemodel MODEL_VERSION=$(MODEL_VERSION)


.PHONY: model-dev
model-dev: test_sta
	mkdir -p $(model_dir)
	rm -f $(model_dir)/*
	cp -v $(model_files) $(model_dir)
	cat modelfiles/lmodel_relxill_devel.dat >> build/lmodel_relxill.dat

	$(eval MODEL_VERSION := $(shell ./test_sta version))
	make compilemodel MODEL_VERSION=$(MODEL_VERSION) DEV=dev

.PHONY: model-nonpublic
model-nonpublic: test_sta
	mkdir -p $(model_dir)
	rm -f $(model_dir)/*
	cp -v $(model_files) $(model_dir)
	cp -v $(unpublished_model_files) $(model_dir)
	cat modelfiles/lmodel_relxill_devel.dat >> build/lmodel_relxill.dat

	$(eval MODEL_VERSION := $(shell ./test_sta version))
	make compilemodel MODEL_VERSION=$(MODEL_VERSION) DEV=nonpublic



.PHONY: valgrind
valgrind:
	make clean
	make CFLAGS="-g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic" test_sta
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./test_sta relconv

valgrind-relxilllp:
	make clean
	make CFLAGS="-g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic" test_sta
	valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all ./test_sta relxilllp

.PHONY: gdb
gdb:
	make clean
	make CFLAGS="-g -ansi -std=c99 -Wall -Wstrict-prototypes -pedantic" test_sta
	gdb --args ./test_sta relline 100

ddd:
	echo "exit" | make gdb 
	ddd test_sta &

gprof:
	make clean
	make CFLAGS="$(CFLAGS) -pg" LDFLAGS="$(LDFLAGS) -pg" test_sta 
	./test_sta relxilllpion 100
	gprof -p test_sta
