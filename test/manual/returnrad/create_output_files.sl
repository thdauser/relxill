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
   set_par("*.logN",_par.logn);     
}
%}}}


define createOutputFiles(_ff, _par){
   putenv("RELXILL_WRITE_FILES=1");

   clean_output_files();
   
%   vmessage("  Creating Output Files: fit-functions: %s ", _ff);

   load_std_model(_ff, _par);
   if (_ff == "relxillBB_NS"){
      set_par("*.shiftTmaxRRad",1.0);
      set_par("*.logN",_par.logn);
   }
   
   %% each evaluation writes the primary (rframe-bbodyEmit, obs-bbodyEmit) and returning (rframe-bbodyRet) spectra
   %% of the zones, and for boost!=0 their reflected spectra, scaled with |boost| (rframe-reflect, obs-reflect, or
   %% rframe-mirror-reflect, obs-mirror-reflect for the mirrored returning spectrum); the total spectrum of a zone is
   %% the sum of its primary and reflected spectrum
   set_par("*.boost",-1,0,-10,10);

   %% 1 creating the "mirror" files
   putenv("RELXILL_RRAD_BBODY_MIRRORED=1");
   () = eval_fun(1,2);

   if (not qualifier_exists("only_mirrored")){
      %% 2 creating the "reflection" files
      putenv("RELXILL_RRAD_BBODY_MIRRORED=0");
      () = eval_fun(1,2);
   }


   putenv("RELXILL_WRITE_FILES=0");
   putenv("RELXILL_RRAD_BBODY_MIRRORED=0");
}
