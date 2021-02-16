require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";


define set_default_dens_cutoff(ff){ %{{{
   if (string_match(ff,"Cp")){
      set_par("*.kTe",30.0);
   }
   
   if (get_params("*.Ecut")[0]!=NULL){
      set_par("*.Ecut",300.0);
   }   
}
%}}}


define check_dens_mod_single(ff,ff_dens){ %{{{
   
%   variable lxi = 2.0;
   
   variable lo,hi;
   variable val0, val1;
   (lo,hi) = log_grid(0.5,500,4000);
   
   if (string_match(ff_dens,"lpD")){
      putenv("RELXILL_NUM_RZONES=1");
   }
   
   fit_fun_default(ff_dens);
   set_par("*.logN",15);
   set_default_dens_cutoff(ff_dens);   
   val0 = eval_fun_keV (lo,hi);
   
   fit_fun_default(ff);
   set_default_dens_cutoff(ff);   
   val1 = eval_fun_keV (lo,hi);

   simple_plot(lo,hi,val0,val1;;__qualifiers());

   putenv("RELXILL_NUM_RZONES");

   variable gn = goodness(val0,val1);
   msg_log += sprintf("    -> %s goodness value: %.3e\n",
	    ff_dens, gn);
   return gn;
}
%}}}



define check_dens_mod_single_shouldChange(ff_dens){ %{{{
   
%   variable lxi = 2.0;
   
   variable lo,hi;
   variable val0, val1;
   (lo,hi) = log_grid(0.5,500,4000);
   
   if (string_match(ff_dens,"lpD")){
      putenv("RELXILL_NUM_RZONES=1");
   }
   
   fit_fun_default(ff_dens);
   set_default_dens_cutoff(ff_dens);
   set_par("*.logN",15);
   val0 = eval_fun_keV (lo,hi);

      
   set_par("*.logN",18);
   val1 = eval_fun_keV (lo,hi);

   val1 *= sum(val0)/sum(val1);
   
   simple_plot(lo,hi,val0,val1;;__qualifiers());

   putenv("RELXILL_NUM_RZONES");

   variable gn = goodness(val0,val1);
   msg_log += sprintf("    -> %s goodness value: %.3e (should fail)\n",
	    ff_dens, gn);
   return gn;
}
%}}}


define runtest(ffs){ 
   
   counter++;
   msg_log += sprintf(" testing HIGH DENSITY MODELS\n");

   variable ff = ["xillver","relxill","relxilllp","xillverCp","relxillCp"];
   variable ff_dens = ["xillverD","relxillD","relxilllpD","xillverDCp","relxillDCp"];

      
   variable ii,n = length(ff);
   
   _for ii(0,n-1,1){
      
      if ( model_exists(ff_dens[ii], ffs)  ){
	 
	 if (check_dens_mod_single_shouldChange(ff_dens[ii];nopl) < goodness_lim*0.1){
	    msg_log += sprintf(" *** error: there seems to be a problem with the HIGH DENSITY MODEL %s\n",ff_dens[ii]);
	    msg_log += sprintf(" ***        -> density parameter does not seem to have an effect\n",ff_dens[ii]); 
	    return EXIT_FAILURE;
	 }
	 
	 if (check_dens_mod_single(ff[ii],ff_dens[ii];nopl) > goodness_lim*0.1){
	    msg_log += sprintf(" *** error: there seems to be a problem with the HIGH DENSITY MODEL %s\n",ff_dens[ii]);
	    msg_log += sprintf(" ***        -> it does not agree with the normal version for logN=15 as it should\n",ff_dens[ii]);
	    return EXIT_FAILURE;
	 }

      }
   }

   

   return EXIT_SUCCESS;
   
}

