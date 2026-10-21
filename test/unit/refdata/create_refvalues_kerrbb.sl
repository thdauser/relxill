% creates reference files with the Xspec model kerrbb for the unit test " disk emission (2b) " in
% tests-bbody-returnrad.cpp: relativistic thin disk (Novikov-Thorne temperature profile) with zero torque at the inner
% edge (eta=0), without spectral hardening (hd=1), self-irradiation (rflag=0), and limb darkening (lflag=0). For each
% spin, the accretion rate is chosen such that the maximal temperature of the disk is 1 keV (Mbh = 10 Msun).
% (energy grid as in the unit tests: 50 log bins from 0.01 to 20 keV)
require("isisscripts");
require("xspec");

variable G = 6.674e-8, c = 2.998e10, sigma_sb = 5.6704e-5, keV_per_K = 8.617333e-8, Msun = 1.989e33;  % cgs
variable Mbh = 10.0, Dbh = 10.0, Tmax = 1.0;

% radius of the ISCO in units of r_g (Bardeen et al. 1972)
define r_isco(a) {
   variable z1 = 1 + (1 - a*a)^(1./3) * ((1 + a)^(1./3) + (1 - a)^(1./3));
   variable z2 = sqrt(3*a*a + z1*z1);
   return 3 + z2 - sign(a) * sqrt((3 - z1) * (3 + z1 + 2*z2));
}

% radial dependence of the Novikov-Thorne flux, F ~ f/x^2, with zero torque at the ISCO and r in units of r_g
% (Page & Thorne 1974, eq. 15n; Li et al. 2005, eq. D13), normalized such that F -> r^-3 for large radii
define flux_nt(r, a) {
   variable x = sqrt(r), x0 = sqrt(r_isco(a));
   variable x1 = 2*cos((acos(a) - PI)/3), x2 = 2*cos((acos(a) + PI)/3), x3 = -2*cos(acos(a)/3);
   variable integ = x - x0 - 1.5*a*log(x/x0)
     - 3*(x1-a)^2/(x1*(x1-x2)*(x1-x3)) * log((x-x1)/(x0-x1))
     - 3*(x2-a)^2/(x2*(x2-x1)*(x2-x3)) * log((x-x2)/(x0-x2))
     - 3*(x3-a)^2/(x3*(x3-x1)*(x3-x2)) * log((x-x3)/(x0-x3));
   return integ / (x^4 * (x^3 - 3*x + 2*a));
}

% maximal temperature [keV] of the disk for the accretion rate mdd [1e18 g/s]: sigma T^4 = 3 G M Mdot / (8 pi r_g^3) F
define tmax_keV(a, mdd) {
   variable M = Mbh * Msun, rg = G * M / c^2;
   variable r = r_isco(a) * [1.0001:4.0:#300001];
   variable F = 3 * G * M * mdd * 1e18 / (8 * PI * rg^3) * flux_nt(r, a);
   return (max(F) / sigma_sb)^0.25 * keV_per_K;
}

variable lo, hi;
(lo, hi) = log_grid(0.01, 20, 50);
fit_fun("kerrbb");
variable a, incl, mdd;
foreach a ([0.0, 0.9, 0.998]) {
   mdd = (Tmax / tmax_keV(a, 1.0))^4;
   foreach incl ([30.0, 60.0]) {
      set_par("kerrbb(1).eta", 0);
      set_par("kerrbb(1).a", a);
      set_par("kerrbb(1).i", incl);
      set_par("kerrbb(1).Mbh", Mbh);
      set_par("kerrbb(1).Mdd", mdd);
      set_par("kerrbb(1).Dbh", Dbh);
      set_par("kerrbb(1).hd", 1);
      set_par("kerrbb(1).rflag", 0);
      set_par("kerrbb(1).lflag", 0);
      set_par("kerrbb(1).norm", 1);
      vmessage(" kerrbb: a=%.3f, incl=%.0f, Mbh=%.1f, Mdd=%.4f -> Tmax=%.6f keV", a, incl, Mbh, mdd, tmax_keV(a, mdd));
      fits_write_binary_table(sprintf("refvalues_kerrbb_a%.3f_i%.0f.fits", a, incl), "REFSPEC",
                              struct{elo=lo, ehi=hi, kerrbb=eval_fun_keV(lo, hi)},
                              struct{SPIN=a, INCL=incl, MBH=Mbh, MDD=mdd, DBH=Dbh, TMAX=Tmax});
   }
}
