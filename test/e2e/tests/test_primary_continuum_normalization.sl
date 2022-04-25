require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";


define check_ns_single(ff){ %{{{

   
   
   fit_fun_default(ff);


   variable lo,hi;
   (lo, hi) = log_grid(0.1,1000,3000);
   
   set_par("*.refl_frac",0.0,0,-10,10);
   variable ener_flux  =  sum(eval_fun_keV(lo,hi)*0.5*(lo+hi));
   
   %% convert to ergs * 10^20
   variable keV2erg = 1.602177e-09;
   variable ener_flux_ergs = ener_flux ;
   
   variable ener_flux_xillver_table = 1e15 / (4.0 * PI * keV2erg * 1e20);

   
   %% Normalization of xillver defined as int_0.1^1000 F_E(ergs)*1e20 = 1e15 / (4*PI)
   %% (see Dauser+2016)   
   
   variable gn = abs(ener_flux_ergs/ener_flux_xillver_table - 1.0);
   msg_log += sprintf("   -> ns model %s  [gn=%.2e]\n",ff,gn);
   

   return gn; 
}
%}}}



define runtest(ffs){ 

   msg_log += sprintf("testing NS-Models\n",counter);

   
   variable ff =      ["xillver", "relxill","xillverNS","relxillNS"];  % %warning, relxilllp does not work correctly, as an additional redshift is applied to the primary spectrum
   
   
   variable ii,n = length(ff);
   
   _for ii(0,n-1,1){

      if ( model_exists(ff[ii], ffs)) {      
	 if (not (check_ns_single(ff[ii]) < goodness_lim)){
	    msg_log += sprintf(" *** error: there seems to be a problem with the PRIMARY CONTINUUM NORMALIZATION in  MODEL %s \n",ff[ii]);
	 return EXIT_FAILURE;
	 }
      } 
	     
   }
   
   return EXIT_SUCCESS;
}
