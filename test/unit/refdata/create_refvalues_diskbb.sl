% creates reference files with the diskbb and ezdiskbb models of Xspec for the unit tests in tests-bbody-returnrad.cpp
% (energy grid as in the unit tests: 50 log bins from 0.01 to 20 keV; diskbb for Tin=1 is in refvalues_diskbb.fits)
require("isisscripts");
require("xspec");
variable lo, hi;
(lo, hi) = log_grid(0.01, 20, 50);
variable T;

fit_fun("diskbb");
set_par("diskbb(1).norm", 1);
foreach T ([0.5, 2.0]) {
   set_par("diskbb(1).Tin", T);
   fits_write_binary_table(sprintf("refvalues_diskbb_Tin%.1f.fits", T), "REFSPEC",
                           struct{elo=lo, ehi=hi, diskbb=eval_fun_keV(lo, hi)});
}

fit_fun("ezdiskbb");
set_par("ezdiskbb(1).norm", 1);
foreach T ([0.5, 1.0, 2.0]) {
   set_par("ezdiskbb(1).T_max", T);
   fits_write_binary_table(sprintf("refvalues_ezdiskbb_Tmax%.1f.fits", T), "REFSPEC",
                           struct{elo=lo, ehi=hi, ezdiskbb=eval_fun_keV(lo, hi)});
}
