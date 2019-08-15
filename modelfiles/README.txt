#########################################
#####  README of the RELXILL model  #####
#########################################

Please read the following instructions carefully in order to ensure a
properly working version of the model. Note that this model is
designed and tested to be used within the X-ray spectral fitting
software ISIS and XSPEC. 

(1) Installing the model:

  - All to execute the compile script by "chmod u+r compile_relxill.csh"

  - Execute "./compile_relxill.csh"

(2) Setting up the model environment:

  - Most importantly the model needs to know where the pre-calculated
    tables are located. By default it assumed that they are in the
    current working directory.

  - The recommended approach is to set the environment variable
    "RELXILL_TABLE_PATH" to the path where the tables are stored. That
    means (for a csh environment) putting a line like the following in
    the .cshrc:
    "setenv RELXILL_TABLE_PATH /home/user/data/relline_tables/"

(3) Loading the model in XSPEC:

  From the directory where the model has been installed, the model can
  simply be loaded inside XSPEC: "lmod relxill .".



(4)  For questions and bug reports, please contact thomas.dauser@fau.de

#########################################
#########################################

RELXILL is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

RELXILL is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
For a copy of the GNU General Public License see
<http://www.gnu.org/licenses/>.

Copyright 2019 Thomas Dauser, Remeis Observatory & ECAP

#########################################
#########################################


