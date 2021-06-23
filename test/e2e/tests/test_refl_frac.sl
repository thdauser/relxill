require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";

   
define ncheck_refl_frac_single(ff){ %{{{
   
   variable val0,val1,valr;
   
   fit_fun_default(ff);

   set_par("*.refl_frac",1,0,-10,10);
   val1 =  eval_fun_keV(lo0,hi0);

   set_par("*.refl_frac",-1,0,-10,10);
   val0 =  eval_fun_keV(lo0,hi0);
   
   return goodness(val1,val0);
}
%}}}

define ncheck_fixReflFrac_working(ff){ %{{{
   
   variable val0,val1,valr;
   
   fit_fun_default(ff);
   
   set_par("*.fixReflFrac",0);
   set_par("*.refl_frac",0.1);
   val1 =  eval_fun_keV(lo0,hi0);

   set_par("*.fixReflFrac",1);
   val0 =  eval_fun_keV(lo0,hi0);
   
   return goodness(val1,val0);
}
%}}}

define ncheck_fixReflFrac_onlyRefl(ff){ %{{{
   
   variable val0,val1,valr;
   
   fit_fun_default(ff);
   
   set_par("*.fixReflFrac",1);
   val1 =  eval_fun_keV(lo0,hi0);

   set_par("*.fixReflFrac",3);
   val0 =  eval_fun_keV(lo0,hi0);
   
   return goodness(val1,val0);
}
%}}}

define check_fixReflFrac(ff){ %{{{
   
   variable val0,val1,valp;
   
   fit_fun_default(ff);

   set_par("*.fixReflFrac",1);
   val1 =  eval_fun_keV(lo0,hi0);

   set_par("*.fixReflFrac",3);
   val0 =  eval_fun_keV(lo0,hi0);

   set_par("*.fixReflFrac",0);
   set_par("*.refl_frac",0);
   valp =  eval_fun_keV(lo0,hi0);

   return goodness(val1,val0+valp);
}
%}}}

define check_refl_frac_single(ff){ %{{{
   
   variable val0,val1,valr;
 
   variable lo0, hi0;
   (lo0, hi0) = log_grid(0.1,10,1000);  %% if we choose something larger than 10, we bet a problem with low numbers from the bbret model
   
   fit_fun_default(ff);
   set_par("*.refl_frac",1.0,0,-10,10);
   val1 =  eval_fun_keV(lo0,hi0);

   set_par("*.refl_frac",0,0,-10,10);
   val0 =  eval_fun_keV(lo0,hi0);

   set_par("*.refl_frac",-1,0,-10,10);
   valr =  eval_fun_keV(lo0,hi0);

   if (qualifier_exists("plot")){
      xlog;%ylog;
      variable ind=where(val1!=0);
      hplot(lo0[ind],hi0[ind],(val0[ind]+valr[ind])/val1[ind]);
      xlog;ylog;
%      hplot(lo0,hi0,val1);
%      ohplot(lo0,hi0,val0);
%      ohplot(lo0,hi0,valr);
%      sleep(100);
   }
   
   return goodness(val1,val0+valr);
}
%}}}

define runtest(ffs){ 
   
   msg_log += sprintf(" testing REFLECTION FRACTION PARAMTER\n");
   
   variable ff = [getAllXillverTypeModels, getAllRelxillTypeModels()];
         
   variable ii,n = length(ff);
   
   variable goodn;
   _for ii(0,n-1,1){
      msg_log += sprintf("   -> reflection fraction in %s\n",ff[ii]);
      goodn = check_refl_frac_single(ff[ii]);
      if (goodn > goodness_lim*5){
	msg_log += sprintf(" *** error: there seems to be a problem with the REFLECTION FRACTION in  MODEL %s (goodness %e)\n",ff[ii],goodn);
	 return EXIT_FAILURE;
      }
      goodn = ncheck_refl_frac_single(ff[ii]);
      if (goodn < goodness_lim){
	 msg_log += sprintf(" *** error: there seems to be a problem with the reflection fraction (no difference between 1 and -1) in  MODEL %s (goodness %e)\n",ff[ii],goodn);
	 return EXIT_FAILURE;
      }
      
      
      if (isLpModel(ff[ii])) {
	 msg_log += sprintf("    + fixReflFrac parameter in %s\n",ff[ii]);
	 if (check_fixReflFrac(ff[ii]) > goodness_lim){
	    msg_log += sprintf(" *** error: there seems to be a problem with the fixReflFrac parameter in  MODEL %s\n",ff[ii]);
	    return EXIT_FAILURE;
	 }
	 if (ncheck_fixReflFrac_working(ff[ii]) < goodness_lim){
	    msg_log += sprintf(" *** error: there seems to be a problem with the fixReflFrac parameter (no difference 0 and 1) in  MODEL %s\n",ff[ii]);
	    return EXIT_FAILURE;
	 }
	 if (ncheck_fixReflFrac_onlyRefl(ff[ii]) < goodness_lim){
	    msg_log += sprintf(" *** error: there seems to be a problem with the fixReflFrac parameter (no difference between 1 and -1) in  MODEL %s\n",ff[ii]);
	    return EXIT_FAILURE;
	 }
      }
   }
   
   return EXIT_SUCCESS;
}
