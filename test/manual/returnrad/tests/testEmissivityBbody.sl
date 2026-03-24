_traceback=1;

require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

variable spin = 0.998;  %% currently only done for this value
variable Tin = 2.0;
variable logxi = 2.;
variable logn = 17.;

define test_primary_emissivity_with_theory(){

   variable fname = getFilenameStruct().rframe.specPri;
   variable all_dat = get_2d_data(fname, spin);
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
   
   variable r = 0.5*(dat.rlo+dat.rhi);
   variable emis_canonical = r^(-3)*(1-1./sqrt(r/kerr_rms(spin)));
   emis_canonical *= ener_flu[-5]/emis_canonical[-5];

   variable ratio = ener_flu/emis_canonical;


   if (qualifier_exists("plot")){
      xlog;ylog;
      plot(r, ener_flu);
      oplot(r, emis_canonical);      
      %plot(r, ener_flu/emis_canonical);
   }

   if (abs(mean(ratio)-1) > 1e-4){
      vmessage(" *** error: testEmissivityBbody() emissivity profile of the inital BBody component does not follow SS73");
      return EXIT_FAILURE;
   }

   vmessage(" \n   --> SUCCESSFUL: emissivity profile of the inital BBody component agrees with SS73 ");
   return EXIT_SUCCESS;	    

}

define test_return_emissivity(){

   variable fname = getFilenameStruct().emis.rrad;

   message(" return emissivity not yet implemented");
   
   return EXIT_SUCCESS;
}

define testEmissivityBbody(){


   load_relxill_model_dev();
   fit_fun("relxillBBxill");
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
