_traceback=1;

require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

variable spin = 0.998;  %% currently only done for this value
variable Tin = 2.0;
variable logxi = 2.;
variable logn = 18.;  %% the only value allowed for relxillBB

define test_primary_emissivity_with_theory(){

   variable fname = getFilenameStruct().rframe.specPri;
   variable all_dat = get_2d_data(fname);
   variable dat = all_dat.dat;
   
   variable jj, m = length(dat.rlo);
   variable ener_flu = Double_Type[m];
   
   _for jj(0,m-1){
      
      % data is in photons/bin
      % we need integrated energy flux / time
      
      % make life easier
      variable elo = all_dat.elo;
      variable ehi = all_dat.ehi;
      variable emn = 0.5*(elo+ehi);
      
      ener_flu[jj] = sum(dat.spec[jj,*]*emn);
      
   }
   
   % canonical emissivity (alpha disk with zero torque at rin, F ~ T^4 ~ r^-3 (1 - sqrt(rin/r))), averaged over the
   % area of each zone, as the primary of a zone is the mean over its area (int F r dr = -1/r + 2/3 sqrt(rin) r^-1.5)
   variable rin = kerr_rms(spin);
   variable r = 0.5*(dat.rlo+dat.rhi);
   variable integ_lo = -1./dat.rlo + 2./3*sqrt(rin)*dat.rlo^(-1.5);
   variable integ_hi = -1./dat.rhi + 2./3*sqrt(rin)*dat.rhi^(-1.5);
   variable emis_canonical = (integ_hi - integ_lo) / (0.5*(dat.rhi^2 - dat.rlo^2));
   emis_canonical *= sum(ener_flu)/sum(emis_canonical);

   variable ratio = ener_flu/emis_canonical;
   variable max_dev = max(abs(ratio/mean(ratio)-1));
   vmessage("   max. deviation of the energy flux of the zones from SS73: %.2e", max_dev);

   if (qualifier_exists("plot")){
      xlog;ylog;
      plot(r, ener_flu);
      oplot(r, emis_canonical);      
      %plot(r, ener_flu/emis_canonical);
   }

   if (max_dev > 1e-3){
      vmessage(" *** error: testEmissivityBbody() emissivity profile of the inital BBody component does not follow SS73");
      return EXIT_FAILURE;
   }

   vmessage(" \n   --> SUCCESSFUL: emissivity profile of the inital BBody component agrees with SS73 ");
   return EXIT_SUCCESS;	    

}

define test_return_emissivity(){

   variable fname = getFilenameStruct().disk.rrad;

   message(" return emissivity not yet implemented");
   
   return EXIT_SUCCESS;
}

define testEmissivityBbody(){


   load_relxill_model_dev();
   fit_fun("relxillBB");
   list_par;
   set_par("*.logxi",logxi);
   set_par("*.logN",logn);
   
   vmessage(" \n Testing testEmissivityBbody() ");
   putenv("RELXILL_WRITE_FILES=1");
   () = eval_fun(1,2);

   variable ntest = 2;
   variable rval = Int_Type[ntest];
   rval[0] = test_primary_emissivity_with_theory();
   rval[1] = test_return_emissivity();

   if (rval[0]!=EXIT_SUCCESS || rval[1]!=EXIT_SUCCESS){
      return EXIT_FAILURE;
   }
   return EXIT_SUCCESS;
}

if (length(__argv)>1){
   () = testEmissivityBbody();
}
