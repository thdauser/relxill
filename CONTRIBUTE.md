# Contributing to the RELXILL model

Contributions to `relxill` are very welcome. Below are a few guidelines 
on how the setup works. General information on the model can be found at
https://www.sternwarte.uni-erlangen.de/research/relxill/

If you are looking for a quick start guide, have a look at the README.md.


####The following topics are covered in this guide:
* general code setup
* the build process
* adding a new model
* test setup 

#### Requirements:
* heasoft (for building the model)
* isis (only for the e2e-tests)

## General Code Setup 

The `relxill` model is written in C and C++. All source files have to be
located in the subdirectory `src/`. Additional files needed to compile the
model (such as the lmodel.dat file), are in `src/modelfiles`. Finally, tests
are put in the folder `test/` (see below for more information on tests).



## The BUILD Process 

The master build process is done via the Makefile in the main
directory. Any compilation or creating of files needed for the
compilation is done in `cmake`, meaning all the information is 
contained in CMakeLists.txt in the main folder, and the respective
subdirectories referenced there.

In order to build all files necessary for the local model, the following steps are performed:
1) build the master `lmodel_relxill.dat` file from the source files 
2) parsing this definition, the `xspec_wrapper_lmodels.cpp/h` files are created, 
   containing the C++ local model definitions
3) then the model files are compiled
4) the model files are collected, and then compiled with `xspec` to build a local model in `build/`

You can test the first three steps, including the compilation, by running `make install`.  

## Adding a new File

In order to add a new source or header file, it has to be added to the `SOURCE_FILES` cmake variable
in the CMakeLists.txt of the current folder. It will then be automatically used for the build process
and if it is in the `src/` directory, it will also be included in the local model and the tar file. Allowed
file endings are `.c .cpp .h`.

## Adding a new Model

A new model is first and foremost defined by the `lmodel.dat` file. How such an entry 
has to be structured is defined 
in Xspec: https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html

1)  First add the model to the `lmodel.dat` file. There are two version of the file contained in 
`src/modelfiles`, a public (`lmodel_relxill_public.dat`) and a development (`lmodel_relxill_devel.dat`) one. 
The difference is that only the public models are included in the official release of the relxill model and 
therefore should only contain stable models with a stable interface.

2) Any new parameters have to be added in `class XPar` and a default value defined  in 
   `class XspecSingleLmodelDefinition` (ModelParams.h)

3) The new model name as to be added as to be added in the `class ModelName` (ModelInfo.h). 
   It is used to uniquely identify the model.
   
4) Additionally, a unique model integer 
   value has to be defined in relmodels.h (for example `#define MOD_TYPE_RELXILL 123`)  and 
   linked to the model, by adding it to the function `int convertModelType(ModelName name)`
   (*this part will soon be removed*).

5) The final model is defined in `class ModelDatabase` (ModelDatabase.h), where the 
   model type (Convolution, Relxill-Type, Line-Model, ...), 
   irradiation type (Lamp Post, Power Law, ...), and
   primary spectral shape (Cutoff-Powerlaw, Nthcomp, Blackbody) are set. 



## Test Setup

The tests are all collected in `test`. The main tests, which are automatically performed, 
are in `test/unit/` and `test/e2e/`. 

From the main directory, calling `make test` will run both, the unit tests and all e2e tests for all
models (stable and development). Additionally, only the model including the stable models can be 
tested as well with `make test-stable`. In this case, first the stable model is built, and then 
the e2e tests are ran on this model. 

#### Unit Tests (test/unit/)
The units tests are done on C or C++ and should test
basic functionality of the code. Those tests are ran first. To implement the tests the Catch2 
framework is used. New tests can simply be added as a new .cpp file and added in CMakeLists.txt.
Those will then be automatically compiled and run.

#### e2e Tests (test/e2e/)
End-to-end tests are based on using the fully compiled local model (in build/), which is then loaded in isis.
Currently, the tests are split in two main functions. All tests are called from the Makefile contained
in this directory. Executing `make test` runs all tests. *Note, these test require ISIS to be installed.*

* `make test-refdata` compares all available reference data to the models which are currently installed. 
Only models for which reference data exists are tested. The script `create_reference_data.sl` can be used 
  to create new reference data.

* `make test-functionality` calls a large isis-script, testing all functions defined in the function
`get_implemented_fitfunctions()`. Note that the models are distinguished in stable and development.
