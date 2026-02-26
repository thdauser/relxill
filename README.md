# relxill - Astrophysics local model for Relativistic Reflection 

## 1. General

“RELXILL” is a X-ray relativistic reflection model developed and maintained by Thomas Dauser, Javier Garica, and collaborators.

The model description and download of the offical version can be found on the homepage at 
http://www.sternwarte.uni-erlangen.de/research/relxill/.

## 2. Scientific Citation 

If you use RELXILL, please cite the relevant papers. 
See: [Name and Citation Policy](NAME_AND_CITATION_POLICY.md)

## 3. Installation
In order to install and build the model from the GIT distribution, simply run "make model". 
It will be installed in the subdirectory "build/". After the installation, a simple test is executed. If this test succeed the model should be correctly working. This installation only works if you have a recent version of Heasoft installed (https://heasarc.gsfc.nasa.gov/lheasoft/).

## 4. Usage
The relxill code is written in C and C++. While its routines can be used directly, its main intention is to be used as *local model* in X-ray data analysis software such as ISIS, Xspec, or Sherpa. The description and meaning of the model parameters can be found in the relxill documentation (https://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/relxill_docu.pdf)
and on the relxill homepage. If the default installation does not work for you, you might need to check how to install local models for your data analysis software.

## 5. Contributing

More information on the code and how to contribute can be found at CONTRIBUTE.md.

## 6. Support
For questions or bug reports, please contact thomas.dauser@fau.de and javier@caltech.edu
