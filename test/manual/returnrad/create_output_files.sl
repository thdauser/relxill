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


%% writes the diagnostic files of the zones (see "Diagnostic output files" in doc/returnrad_disk_emission.md) for the
%% reflection only (boost=-1), once with the mirrored returning radiation (rframe-mirror-reflect, obs-mirror-reflect) and
%% once with the xillver reflection (rframe-reflect, obs-reflect); both write the primary (rframe-bbodyEmit,
%% obs-bbodyEmit) and the returning radiation (rframe-bbodyRet), which are the same in both cases
%% qualifier order: order of the returning radiation (RELXILL_RRAD_BBODY_ORDER, default 1)
define createOutputFiles(_ff, _par){
   putenv("RELXILL_WRITE_FILES=1");
   putenv(sprintf("RELXILL_RRAD_BBODY_ORDER=%d", qualifier("order",1)));

   clean_output_files();

   load_std_model(_ff, _par);
   set_par("*.boost",-1,0,-10,10);

   putenv("RELXILL_RRAD_BBODY_MIRRORED=1");
   () = eval_fun(1,2);

   putenv("RELXILL_RRAD_BBODY_MIRRORED=0");
   () = eval_fun(1,2);

   putenv("RELXILL_WRITE_FILES=0");
   putenv("RELXILL_RRAD_BBODY_ORDER=1");
}
