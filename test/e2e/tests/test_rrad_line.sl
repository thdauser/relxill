require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";

define single_comparison(ff, ff_ref){
   

   variable lo, hi;
   (lo,hi) = log_grid(0.2,10,600);
   

   fit_fun(ff_ref);
   set_par("*.lineE",6.4);
   variable val_ref = eval_fun_keV(lo,hi);

   fit_fun(ff);
   set_par("*.lineE",6.4);
   set_par("*return_rad",qualifier("return_radiation"));
   variable val0 = eval_fun_keV(lo,hi);
   
   variable i_en = where(2.0<lo<6.0);
   variable gn = goodness(val0[i_en],val_ref[i_en]);
   msg_log += sprintf("    -> %s compared with %s goodness value: %.3e (rrad=%i)\n", ff, ff_ref, gn,qualifier("return_radiation"));
   
   if ( qualifier_exists("pl")){
      hplot(lo,hi,val0);
      ohplot(lo,hi,val_ref);
   }

   return gn;
}

define runtest(ffs){

   variable status = EXIT_SUCCESS;

   variable ff = "rellinelpRet";
   variable ff_ref = "relline_lp";

   if (not model_exists(ff,ffs)){
      message(" ... skipping ");
      return status;
   }
   

   
   if ( single_comparison(ff, ff_ref; return_radiation=0) > goodness_lim){
      msg_log += sprintf(" *** error: matching rrad line profile to standard profile\n");
      status = EXIT_FAILURE;
   }
   if ( single_comparison(ff, ff_ref; return_radiation=1) < goodness_lim){
      msg_log += sprintf(" *** error: rrad line profile does not differ from standard profile\n");
      status = EXIT_FAILURE;
   }
   
   return status;
}
