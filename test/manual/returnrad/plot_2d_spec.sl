require("isisscripts");


variable fname = "cmake-build-debug/debug-testrr-spec-diskbb.fits";
variable fname_spec = "cmake-build-debug/testrr-spec-diskbb.fits";

define plot_dat_2d(fname, fname_spec){
   
   
   variable dat = fits_read_table (fname);
   
   variable ii, n = length(dat.rlo);
   
   variable pl = xfig_plot_new();
   pl.world(0.05,20,1e-3,100;loglog);
   variable sumspec = Double_Type[n];
   variable areafac = Double_Type[n];
   _for ii(0, n-1){

      areafac[ii] = (dat.rhi[ii]^2-dat.rlo[ii]^2)/(dat.rhi[n-1]^2 - dat.rlo[0]^2);
      
%      dat.spec[ii,*] *= areafac[ii];
      
      pl.plot(dat.ener[ii,*], dat.spec[ii,*]);


      sumspec += dat.spec[ii,*];
   }
   
   
   print(sum(areafac));
   
   pl.plot(dat.ener[0,*], sumspec; width=3);
   

   variable spec = fits_read_table(fname_spec); %% ener-array is actually bin_lo
   pl.plot(spec.bin_lo, spec.flux / sum(spec.flux) * sum(sumspec);
	  width=3, color="blue");

   
   fit_fun("diskbb");
   set_par("*.Tin",1.0);
   variable lo = spec.bin_lo;
   variable hi = spec.bin_hi;
   variable val = eval_fun_keV(lo,hi);
   
   pl.plot(lo, val/sum(val)*sum(sumspec); width=3, color="red");

   
   return pl;
}


variable pl = plot_dat_2d(fname,fname_spec);


pl.render("plot_diskbb_comparison.pdf");