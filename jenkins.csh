#!/bin/csh

if (-d /data/system/software ) then
    setenv SOFTDIR /data/system/software
#    source /data/system/software/softwarescript.csh                                                                     
     /data/system/software/softwarescript_Xray.csh                                                                
else
    setenv SOFTDIR /home/${USER}/software
endif

if (-d $SOFTDIR/softwarescript_Xray.csh ) then
    source $SOFTDIR/softwarescript_Xray.csh
endif

echo DONE

#setenv PATH "$PATH:/data/system/software/isis/current/x86_64-libc2.19/bin"
#setenv PATH "$PATH:/data/system/software/slang/x86_64-libc2.15/bin"

export RELXILL_TABLE_PATH="/userdata/data/dauser/relline_tables"

make 
./test_sta

make model

cd test/ && ./check_model_functions.sl
