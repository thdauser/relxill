_traceback=1;

require("isisscripts");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

define testEmissivityBbody(){

   variable spin = 0.998;  %% currently only done for this value
   variable Tin = 2.0;
   
   fit_fun("relxillBBxill");
   
   vmessage(" \n Testing testEmissivityBbody() ");
   putenv("RELXILL_WRITE_FILES=1");
   () = eval_fun(1,2);
   
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

if (length(__argv)>1){
   () = testXillverPrimaryNormalization();
}
