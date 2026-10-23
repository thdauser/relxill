#!/usr/bin/env isis-script

%% observed spectra of the zones of relxillBB_NS with the mirrored returning radiation (RELXILL_RRAD_BBODY_MIRRORED=1,
%% reflection only): direct emission, mirrored 1st order, and mirrored 2nd order of the returning radiation
%% (plots/relxillBB_NS_perzone_mirror_orders.pdf)
%% The model returns the sum of all orders up to RELXILL_RRAD_BBODY_ORDER, so the 2nd order alone is the difference of
%% the runs with order 2 and order 1. See "Diagnostic output files" in doc/returnrad_disk_emission.md.

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

   %% the direct emission does not depend on the order (written by both runs)
   variable d_direct = get_2d_data(getFilenameStruct().fobs.primary);
   variable d_order1 = get_2d_data(get_fname_mirror_order(1));
   variable d_order12 = get_2d_data(get_fname_mirror_order(2));

   %% 2nd order alone
   variable dat_order2 = @d_order12.dat;
   dat_order2.spec = d_order12.dat.spec - d_order1.dat.spec;
   variable d_order2 = struct{elo=d_order12.elo, ehi=d_order12.ehi, sumspec=d_order12.sumspec - d_order1.sumspec,
			      dat=dat_order2};

   variable title = "Observed, mirror: 1st order (solid), 2nd order (dashed), direct (dotted); "
		   + "sum of the zones (black), direct + 1st + 2nd order (red)";
   variable pl = get_2d_plot_data([d_direct, d_order1, d_order2], spin; title=title, total=[0,1,2], lines=[2,0,1]);
   pl.render(dir_plots+ff_plots+"_perzone_mirror_orders.pdf");
}

plot_perzone_mirror_orders();
