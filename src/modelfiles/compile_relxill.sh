#!/bin/bash

ln -sf lmodel_relxill.dat lmodel.dat

sed -i "s,#define RELXILL_TABLE_PATH.*,#define RELXILL_TABLE_PATH "'"'`pwd`'"'"," common.h

#  "For Mac OSX systems, you might have to use this"
# sed -i.ori "s,#define RELXILL_TABLE_PATH.*,#define RELXILL_TABLE_PATH "'"'`pwd`'"'"," common.h

echo "initpackage relxill lmodel_relxill.dat `pwd` \nexit" | xspec
echo "lmod relxill . \nexit" | xspec

rm -f *~ *.o *FunctionMap.* lpack_* *.mod Makefile
