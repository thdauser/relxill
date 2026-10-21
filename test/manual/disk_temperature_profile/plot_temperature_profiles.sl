#!/usr/bin/env isis-script
% Plot of the temperature profiles of the disk emission of the relxillBB models (see doc/returnrad_disk_emission.md):
% the alpha disk with zero torque at the inner edge (default, RELXILL_RRAD_BBODY_TPROFILE=1) and the Novikov-Thorne
% profile (RELXILL_RRAD_BBODY_TPROFILE=3) for several spins, both normalized to their maximal temperature (= kTbb),
% with the inner edge of the disk at the ISCO.
%
% The profiles are calculated with the same equations as in src/Relphysics.cpp (disk_Tprofile_alpha() and
% disk_Tprofile_novikov_thorne()).
%
% Usage:  isis-script plot_temperature_profiles.sl
%         (writes disk_temperature_profiles.pdf in this directory and doc/images/disk_temperature_profiles.png)

require("isisscripts");

variable spins = [0.0, 0.9, 0.998];

% radius of the ISCO in units of r_g (Bardeen et al. 1972)
define r_isco(a) {
   variable z1 = 1 + (1 - a*a)^(1./3) * ((1 + a)^(1./3) + (1 - a)^(1./3));
   variable z2 = sqrt(3*a*a + z1*z1);
   return 3 + z2 - sign(a) * sqrt((3 - z1) * (3 + z1 + 2*z2));
}

% alpha disk with zero torque at rin, normalized to T = 1 at Rmax = 1.5^(4/5) rin (as disk_Tprofile_alpha())
define tprofile_alpha(r, rin) {
   variable rmax = 1.5^0.8 * rin;
   return (r/rin)^(-0.75) * (1 - sqrt(rin/r))^0.25 / ((rmax/rin)^(-0.75) * (1 - sqrt(rin/rmax))^0.25);
}

% radial dependence of the Novikov-Thorne flux with zero torque at rin, F ~ f/x^2 (Page & Thorne 1974, eq. 15n;
% Li et al. 2005, eq. D13), with r in units of r_g (as disk_flux_novikov_thorne())
define flux_nt(r, rin, a) {
   variable x = sqrt(r), x0 = sqrt(rin);
   variable x1 = 2*cos((acos(a) - PI)/3), x2 = 2*cos((acos(a) + PI)/3), x3 = -2*cos(acos(a)/3);
   variable integ = x - x0 - 1.5*a*log(x/x0)
     - 3*(x1-a)^2/(x1*(x1-x2)*(x1-x3)) * log((x-x1)/(x0-x1))
     - 3*(x2-a)^2/(x2*(x2-x1)*(x2-x3)) * log((x-x2)/(x0-x2))
     - 3*(x3-a)^2/(x3*(x3-x1)*(x3-x2)) * log((x-x3)/(x0-x3));
   integ *= (integ > 0);   % only rounding can make it negative
   return integ / (x^4 * (x^3 - 3*x + 2*a));
}

% Novikov-Thorne profile with zero torque at rin, normalized to its maximal temperature
define tprofile_nt(r, rin, a) {
   variable rfine = rin * [1.000001:4.0:#200001];
   return (flux_nt(r, rin, a) / max(flux_nt(rfine, rin, a)))^0.25;
}


%%% plot %%%
xfig_new_color("ink2", 0x52514e);
variable cols = ["nt_a0", "nt_a09", "nt_a0998"];
xfig_new_color(cols[0], 0x2a78d6);
xfig_new_color(cols[1], 0xeb6834);
xfig_new_color(cols[2], 0x1baf7a);

variable x = 10^[log10(1.000001):3.0:#3000];   % r / Rin
variable t_alpha = tprofile_alpha(x, 1.0);       % independent of the spin in units of r / Rin

variable pl = xfig_plot_new(14, 8);
pl.world(1, 1000, 0, 1.08; xlog);
pl.ylabel("$T / T_\mathrm{max}$ \ ($T_\mathrm{max} = kT_\mathrm{bb}$)"R);
pl.title("Disk temperature profiles, inner edge at the ISCO");

variable pl_ratio = xfig_plot_new(14, 5);
pl_ratio.world(1, 1000, 0.82, 1.28; xlog);
pl_ratio.xlabel("$r / R_\mathrm{in}$"R);
pl_ratio.ylabel("$T_\mathrm{NT} / T_\mathrm{alpha}$"R);
pl_ratio.plot([1, 1000], [1, 1]; color="ink2", width=1);
pl_ratio.xylabel(1.15, 1.25, "Novikov-Thorne relative to the alpha disk (same $T_\mathrm{max}$)"R, -0.5, 0;
                 color="ink2", size="\footnotesize"R);

variable ii, imax, labels = String_Type[length(spins) + 1];
pl.plot(x, t_alpha; color="ink2", width=3, line=1);
imax = where_max(t_alpha)[0];
pl.plot(x[imax], 1.0; sym="circle", fill=20, size=0.4, color="ink2");
labels[0] = sprintf("alpha disk (any spin), $R_\mathrm{max} = %.2f\,R_\mathrm{in}$"R, x[imax]);

_for ii (0, length(spins) - 1) {
   variable a = spins[ii], rin = r_isco(a);
   variable t_nt = tprofile_nt(x * rin, rin, a);
   pl.plot(x, t_nt; color=cols[ii], width=3);
   imax = where_max(t_nt)[0];
   pl.plot(x[imax], 1.0; sym="circle", fill=20, size=0.4, color=cols[ii]);
   labels[ii + 1] = sprintf("Novikov-Thorne, $a = %g$, $R_\mathrm{max} = %.2f\,R_\mathrm{in}$"R, a, x[imax]);

   variable m = where(x > 1.02);
   variable ratio = t_nt[m] / t_alpha[m];
   pl_ratio.plot(x[m], ratio; color=cols[ii], width=3);
   pl_ratio.xylabel(x[m][-1], ratio[-1], sprintf(" $a = %g$", a), -0.5, 0; size="\footnotesize"R);
}
pl.add_object(xfig_new_legend(labels, ["ink2", cols], [1, 0, 0, 0], 3, 0.8; labelsize="\footnotesize"R),
              0.97, 0.95, 0.5, 0.5; world0);

variable script_dir = path_dirname(__FILE__);
variable fname_pdf = path_concat(script_dir, "disk_temperature_profiles.pdf");
variable fname_png = path_concat(script_dir, "../../../doc/images/disk_temperature_profiles.png");
xfig_multiplot(pl, pl_ratio; cols=1).render(fname_pdf);
() = system(sprintf("convert -density 200 -background white -alpha remove %s %s", fname_pdf, fname_png));
vmessage(" plot written to %s and %s", fname_pdf, fname_png);
