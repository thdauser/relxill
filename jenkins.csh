#!/bin/csh

if (-d /data/system/software ) then
    setenv SOFTDIR /data/system/software
    setenv RELXILL_TABLE_PATH "/userdata/data/dauser/relline_tables"
else if (-d /home/thomas/software) then
    setenv SOFTDIR /home/thomas/software
    setenv RELXILL_TABLE_PATH "/home/thomas/data/relline_tables"
endif

if (-e $SOFTDIR/softwarescript_Xray.csh ) then

    if ($HOST == "cepheus") then
	setenv PATH /data/system/software/isis/current/x86_64-libc2.23/bin:$PATH
	echo "Using isis in /data/system/software/isis/current/x86_64-libc2.23/bin"
    else if (-e $SOFTDIR/softwarescript.csh ) then
	source $SOFTDIR/softwarescript.csh
    endif
    source $SOFTDIR/softwarescript_Xray.csh
    echo " Using scripts $SOFTDIR/softwarescript_Xray.csh "
endif

make clean

make
./test_sta > test_output.log

make model
test/check_model_functions.sl  >> test_output.log
