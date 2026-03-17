require("isisscripts");
require("scripts/subs_returnrad.sl");

variable dir = "build/";

variable fname = "cmake-build-debug/testrr-spec-diskbb.fits";
variable fd = fits_read_table(fname);

variable ref = fits_read_table("refdata/refvalues_diskbb.fits");

variable elo = fd.bin_lo;
variable ehi = fd.bin_hi;
variable emn = 0.5*(elo+ehi);
variable nener = length(elo);

variable pl = xfig_plot_new();
variable plr = xfig_plot_new();

pl.world(0.01,20,1e-2,50;loglog);
plr.world(0.01,20,0.3,1.3;xlog);

fd.flux *= emn/(ehi-elo);
pl.plot(emn,fd.flux);

fit_fun("diskbb");
set_par("diskbb(1).Tin",1.0);

variable val = eval_fun_keV(elo,ehi)/(ehi-elo)*emn;

val *= fd.flux[nener/2]/val[nener/2];

pl.plot(emn,val;color="blue");

ref.spec *= 2*PI;
pl.plot(0.5*(ref.elo+ref.ehi),ref.spec;color="red");



plr.plot(emn,fd.flux/val;color="blue");
plr.plot(emn,fd.flux/ref.spec;color="red");

xfig_multiplot(pl,plr).render("test_diskbb.pdf");