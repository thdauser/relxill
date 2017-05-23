#!/bin/tcsh

ln -sf lmodel_relxill.dat lmodel.dat

sed -i "s,#define RELXILL_TABLE_PATH.*,#define RELXILL_TABLE_PATH "'"'`pwd`'"'"," relbase.h

echo "initpackage relxill lmodel_relxill.dat `pwd` \nquit\ny" | xspec
echo "lmod relxill .\nquit\ny" | xspec

rm -f *~ *.o *FunctionMap.* lpack_* *.mod Makefile pkgIndex.tcl
