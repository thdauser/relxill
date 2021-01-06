## relxill - Astrophysics local model for Relativistic Reflection 

Copyright 2021 Thomas Dauser, Remeis Observatory & ECAP

### 1. General
The model description and download of released version can be found on the homepage at 
http://www.sternwarte.uni-erlangen.de/research/relxill/  -- In order to obtain a stable version this
is the best place to get the model from.

### 2. Installation
In order to install and build the model from the GIT distribution, simply run "make model". 
It will be installed in the subdirectory "build/". After the installation, a simple test is executed. If
this test succeed the model should be correctly working. This installation only works if you have a recent
version of Heasoft installed (https://heasarc.gsfc.nasa.gov/lheasoft/).

### 3. Usage
The relxill code is written in C and C++. While its routines can be used directly, its main intention
is to be used as *local model* in X-ray data analysis software such as ISIS, Xspec, or Sherpa. The 
description and meaning of the model parameters can be found on the relxill homepage (see above). If the
default installation does not work for you, you might need to check how to install local models for your
data analysis software.

### 4. Contributing

More information on the code and how to contribute can be found at CONTRIBUTE.md.

### 5. Support
For questions or bug reports, please contact thomas.dauser@fau.de and javier@caltech.edu
