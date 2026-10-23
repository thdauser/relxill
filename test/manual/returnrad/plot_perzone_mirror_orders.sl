#!/usr/bin/env isis-script

%% observed spectra of the zones of relxillBB_NS with the mirrored returning radiation (RELXILL_RRAD_BBODY_MIRRORED=1,
%% reflection only), for each zone and integrated over the zones: direct emission, mirrored 1st order, and mirrored
%% 1st + 2nd order of the returning radiation (plots/relxillBB_NS_perzone_mirror_orders.pdf)
%% The model returns the sum of all orders up to RELXILL_RRAD_BBODY_ORDER, i.e., the run with order 2 contains the 1st
%% and the 2nd order. See "Diagnostic output files" in doc/returnrad_disk_emission.md.

_traceback=1;

require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

require("create_output_files.sl");

variable ff_plots = "relxillBB_NS";

variable spin = 0.998;
variable Tin = 2.0;
variable logxi = 2.0;  %% the mirror does not depend on logxi and logN
variable logn = 18.0;

%% name of the copy of the observed mirrored spectra of the zones of the run with the given order
define get_fname_mirror_order(order){
   return sprintf("./__mirror-order%d-obs-mirror-reflect.fits", order);
}

%% writes the files with the mirrored returning radiation up to the given order, and keeps a copy of the observed
%% mirrored spectra of the zones (the next run overwrites the files)
define create_mirror_files(order){
   putenv("RELXILL_WRITE_FILES=1");
   putenv(sprintf("RELXILL_RRAD_BBODY_ORDER=%d", order));
   putenv("RELXILL_RRAD_BBODY_MIRRORED=1");
   clean_output_files();

   load_std_model(ff_plots, get_param_struct(spin, Tin, logxi, logn));
   set_par("*.boost",-1,0,-10,10);
   () = eval_fun(1,2);

   () = system(sprintf("cp %s %s", getFilenameStruct().fobs.mirrorRefl, get_fname_mirror_order(order)));

   putenv("RELXILL_WRITE_FILES=0");
   putenv("RELXILL_RRAD_BBODY_ORDER=1");
   putenv("RELXILL_RRAD_BBODY_MIRRORED=0");
}

define plot_perzone_mirror_orders(){
   create_mirror_files(1);
   create_mirror_files(2);

   %% the direct emission does not depend on the order (written by both runs); the run with order 2 contains the
   %% 1st and the 2nd order
   variable fdirect = getFilenameStruct().fobs.primary;
   variable title = "Observed, mirror: direct (dotted), 1st order (dashed), 1st + 2nd order (solid); "
		   + "zones (colors) and integrated over the zones (black)";
   variable pl = get_2d_plot({fdirect, get_fname_mirror_order(1), get_fname_mirror_order(2)}, spin;
			     title=title, lines=[2,1,0]);
   pl.render(dir_plots+ff_plots+"_perzone_mirror_orders.pdf");
}

plot_perzone_mirror_orders();
