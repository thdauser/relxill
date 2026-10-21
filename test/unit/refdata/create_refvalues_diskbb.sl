% creates reference files with the diskbb model of Xspec for the unit tests in tests-bbody-returnrad.cpp
% (energy grid as in the unit tests: 50 log bins from 0.01 to 20 keV; Tin=1 is in refvalues_diskbb.fits)
require("isisscripts");
require("xspec");
variable lo, hi;
(lo, hi) = log_grid(0.01, 20, 50);
fit_fun("diskbb");
set_par("diskbb(1).norm", 1);
variable Tin;
foreach Tin ([0.5, 2.0]) {
   set_par("diskbb(1).Tin", Tin);
   fits_write_binary_table(sprintf("refvalues_diskbb_Tin%.1f.fits", Tin), "REFSPEC",
                           struct{elo=lo, ehi=hi, diskbb=eval_fun_keV(lo, hi)});
}
