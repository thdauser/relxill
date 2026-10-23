#!/usr/bin/env isis-script

%% diagnose plots of the spectra of the zones of relxillBB (first order of the returning radiation, reflection only):
%%  - rest frame: direct, returning, and reflected spectrum of the zones, normalized per zone by N_xill of the returning
%%    spectrum (plots/relxillBB_perzone_rframe.pdf)
%%  - observed: direct, mirrored, and reflected spectrum of the zones, and these spectra integrated over the zones
%%    (plots/relxillBB_perzone_observed.pdf)
%% see "Diagnostic output files" in doc/returnrad_disk_emission.md

_traceback=1;

require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

require("create_output_files.sl");

variable ff_plots = "relxillBB";
variable prefix = ff_plots+"_";

%% parameters (the only ones of the xillverRTR table v0.3)
variable spin = 0.998;
variable Tin = 2.0;
variable logxi = 2.0;
variable logn = 18.0;


define plot_diagnose_rframe(){
   variable files = getFilenameStruct();
   variable title = "Rest frame, per zone normalized to $N_\mathrm{xill}$ of the return: "R
		   + "reflected (solid), return (dashed), direct (dotted)";

   %% the direct emission is about 1/(returning fraction) ~ 25 times above the returning spectrum
   variable pl = get_2d_plot({files.rframe.specPri, files.rframe.specRet, files.rframe.xillRefl}, spin;
			     rframe, title=title, ymax=1e5, lines=[2,1,0], nzones=6);
   pl.render(dir_plots+prefix+"perzone_rframe.pdf");
}

define plot_diagnose_observed(){
   variable files = getFilenameStruct();
   variable title = "Observed: direct (dotted), mirrored (dashed), reflected (solid); "
		   + "zones (colors) and integrated over the zones (black)";

   variable pl = get_2d_plot({files.fobs.primary, files.fobs.mirrorRefl, files.fobs.reflect}, spin;
			     title=title, lines=[2,1,0], nzones=6);
   pl.render(dir_plots+prefix+"perzone_observed.pdf");
}


createOutputFiles(ff_plots, get_param_struct(spin, Tin, logxi, logn));
plot_diagnose_rframe();
plot_diagnose_observed();
