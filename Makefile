# -*- mode: Make -*-

BUILD_DIR = "cmake-build"
MODEL_DIR = "model"
MODEL_BUILD_DIR = "model-build"
BIN_DIR = "bin"

COMPILE_SCRIPT = "compile_relxill.sh"

MODEL_VERSION = "undef"
LMODELDAT_TARGET = "lmodeldat"

TARFILE = "relxill.tgz"

.PHONY: model-stable, model-dev model-build-dev, model-build-target, model-compile, model-tarball
model:
	make model-build-dev

model-stable:
	make model-build-stable



.PHONY: model-build-stable
model-build-stable: export RELXILL_STABLE=
model-build-stable:
	make model-build-target DEV=

model-build-dev:
	make model-build-target DEV=dev

model-build-target:
	make clean
	make install-source-files

	$(eval MODEL_VERSION := `$(BIN_DIR)/test_sta version`)
	$(eval MODEL_TAR_NAME := relxill_model_v$(MODEL_VERSION)$(DEV).tgz)

	make model-tarball TARFILE=$(MODEL_TAR_NAME)
	make model-compile TARFILE=$(MODEL_TAR_NAME)

	@echo "\n  --> Built model  *** $(MODEL_TAR_NAME) *** \n"



model-tarball:
	cd $(MODEL_DIR)/ && tar cfvz $(TARFILE) *.c *.h lmodel_relxill.dat $(COMPILE_SCRIPT) -C ../ README.txt LICENSE
	cp -v $(MODEL_DIR)/$(TARFILE) .


model-compile:
	mkdir -p $(MODEL_BUILD_DIR)
	rm -f $(MODEL_BUILD_DIR)/*

	mv $(TARFILE) $(MODEL_BUILD_DIR)/
	cd $(MODEL_BUILD_DIR) && tar xfvz $(TARFILE)
	cd $(MODEL_BUILD_DIR) && chmod a+x $(COMPILE_SCRIPT)
	cd $(MODEL_BUILD_DIR) && ./$(COMPILE_SCRIPT)
	cd $(MODEL_BUILD_DIR) && echo 'require("xspec"); load_xspec_local_models("./librelxill.so"); fit_fun("relxill"); () = eval_fun(1,2); exit; ' | isis -v



# using CMakeLists.txt Definitions to build and install the model files
.PHONY: install, install-stable, install-source-files, build

install-stable: export RELXILL_STABLE=
install-stable:
	make install-source-files

install:
	make install-source-files

build:
	mkdir -p $(BUILD_DIR)
	cd $(BUILD_DIR) && cmake ../
	cd $(BUILD_DIR) && cmake --build .


install-source-files:
	make build
	cd $(BUILD_DIR) && cmake --install .


.PHONY: test
test:
	make test-unit
	make test-e2e


test-unit:
	make install
	bin/run-tests

test-e2e:
	cd test/e2e && make test

test-stable:
	make model-stable
	cd test/e2e && make test-stable


clean:
	rm -f *~ gmon.out test*.dat *.log
	rm -rf $(MODEL_DIR)
	rm -rf $(MODEL_BUILD_DIR)
	rm -rf $(BIN_DIR)
	rm -rf $(BUILD_DIR)
	rm -f relxill_model_v*.tgz
	rm -f debug-*fits testrr-*.fits
	rm -rf cmake-*


dist-clean:
	make clean
	rm -rf cmake-*
