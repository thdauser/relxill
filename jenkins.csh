#!/bin/csh

if (-d /data/system/software ) then
    setenv SOFTDIR /data/system/software
else if (-d /home/${USER}/software)
    setenv SOFTDIR /home/${USER}/software
endif

if (-e $SOFTDIR/softwarescript_Xray.csh ) then
    if (-e $SOFTDIR/softwarescript.csh ) then
	source $SOFTDIR/softwarescript.csh
    endif
    source $SOFTDIR/softwarescript_Xray.csh
    echo " Using scripts $SOFTDIR/softwarescript_Xray.csh "
endif

setenv RELXILL_TABLE_PATH "/userdata/data/dauser/relline_tables"

make 
 ./test_sta

make model

cd test/ && ./check_model_functions.sl
