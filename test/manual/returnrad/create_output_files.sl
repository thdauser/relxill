#!/usr/bin/env isis-script

_traceback=1;

require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

load_xspec_local_models("build/");

define load_std_model(_ff,_par){
%{{{

   fit_fun(_ff);

   set_par("*.a",_par.spin);
   set_par("*.kTbb", _par.Tin);
   set_par("*.Incl",_par.incl);
   set_par("*.logxi",_par.logxi);     
}
%}}}


define createOutputFiles(_ff, _par){
   putenv("RELXILL_WRITE_FILES=1");

   clean_output_files();
   
%   vmessage("  Creating Output Files: fit-functions: %s ", _ff);

   load_std_model(_ff, _par);
   if (_ff == "relxillBB"){
      set_par("*.shiftTmaxRRad",1.0);
      set_par("*.logN",_par.logn);
   }
   
   %% creating the "mirror" files
   putenv("RELXILL_RRAD_BBODY_MIRRORED=1");

   set_par("*.boost",1.0);
   () = eval_fun(1,2);

   set_par("*.boost",0);
   () = eval_fun(1,2);

   set_par("*.boost",-1.0,0,-10,10);
   () = eval_fun(1,2);
   

   if (not qualifier_exists("only_mirrored")){
      
      %%% 2 creating the "reflection" files
      putenv("RELXILL_RRAD_BBODY_MIRRORED=0");
      variable files = getFilenameStruct();
      
      %%% 2a combined
      set_par("*.boost",1);
      () = eval_fun(1,2);
      
      %%% 2b reflected   
      set_par("*.boost",-1.0,0,-10,10);
      () = eval_fun(1,2);
      variable cmd = sprintf("mv %s %s",files.rframe.xillCombined,files.rframe.xillRefl);
      %   message(cmd);
      () = system(cmd);
      
      %%% 2c primary   
      set_par("*.boost",0.0,0,-10,10);
      () = eval_fun(1,2);
      cmd = sprintf("mv %s %s",files.rframe.xillCombined,files.rframe.xillPrim);
      % message(cmd);
      () = system(cmd);
      
      
      set_par("*.boost",0.0,0,-10,10);
      () = eval_fun(1,2);
   }


   putenv("RELXILL_WRITE_FILES=0");
   putenv("RELXILL_RRAD_BBODY_MIRRORED=0");
}
