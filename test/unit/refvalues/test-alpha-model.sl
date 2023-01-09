% -*- mode: slang; mode: fold -*-

require("isisscripts");
require("physical_disk_irradiation");


variable a = 0.998;
variable h = 6.0;
variable rin = kerr_rms(a);  
variable r_lxi_max = rin*(11./9.)^2;  %% radius where lxi==max(lxi)
variable rout = 400.0;
variable gamma = 2.0;

variable mass_msolar = 20.;
variable flux0 = 1e-10;
variable D_kpc = 10.0;

variable lobs = convertFobs2Lobs(flux0, D_kpc);
variable lobs_default = 1e38;

variable setup = struct{
   a = a,
   h = h,
   Rin = rin,
   Rout = rout,
   gamma = gamma,
   msolar = mass_msolar,
   lobs = lobs
};



define get_phys_emis_profile(setup){ %{{{
%%%

   variable physical_irradiation = calc_physical_irradiation(setup, setup.lobs, setup.msolar);


   variable flux_r_lxi_max = calcFluxRin(physical_irradiation.Lsrc, setup.h, r_lxi_max, setup.msolar);
   vmessage(" - L_source = %.4e", physical_irradiation.Lsrc);
   vmessage(" Newtowning Flux at r_lxi_max=%.4f = %.4e",  r_lxi_max, flux_r_lxi_max);

   
   variable emis_profile, emis_profile_newton;
   (emis_profile, emis_profile_newton) = getEmisPhysnorm_NEW(physical_irradiation.model_params, physical_irradiation);
     
   return emis_profile, emis_profile_newton;
}
%}}}


variable emis_profile, emis_profile_newton;
(emis_profile, emis_profile_newton) = get_phys_emis_profile(setup);

xlog;ylog;

variable ind_lxi = wherefirst(emis_profile.rad >= r_lxi_max);
variable emis_lxi_max = emis_profile.emis[ind_lxi];
variable emis_newton_lxi_max = emis_profile_newton[ind_lxi];

plot(emis_profile.rad, emis_profile.emis);
oplot(emis_profile.rad, emis_profile_newton);

variable ind_lxi_hi = wherefirst(emis_profile.rad >= r_lxi_max)-1;
variable emis_lxi_max_hi = emis_profile.emis[ind_lxi_hi];

vmessage("flux at r_lxi_max=[%.3e-%.3e], for radial bin r=[%.3f-%.3f]",
	 emis_lxi_max, emis_lxi_max_hi,
	 emis_profile.rad[ind_lxi],
	 emis_profile.rad[ind_lxi_hi]);


