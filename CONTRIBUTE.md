# Contributing to the RELXILL model

Contributions to `relxill` are very welcome. Below are a few guidelines 
on how the setup works. General information on the model can be found at
https://www.sternwarte.uni-erlangen.de/research/relxill/

If you are looking for a quick start guide, have a look at the README.


The following topics are covered in this guide:
* general code setup
* the build process
* adding a new model
* test setup 


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

2) Any new parameters have to be added and initialized

3) The model type has to be set (will be removed)

4) The final model is defined in ** bla **

## Test Setup