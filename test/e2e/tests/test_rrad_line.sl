require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";

define single_comparison_rrad(ff){
   

   variable lo, hi;
   (lo,hi) = log_grid(0.2,10,600);
   

   fit_fun(ff);
   set_par("*.lineE",6.4);
   set_par("*.switch_returnrad", 0);
   variable val_ref = eval_fun_keV(lo,hi);

   set_par("*switch_returnrad",1);
   variable val0 = eval_fun_keV(lo,hi);
   
   variable i_en = where(2.0<lo<6.0);
   variable gn = goodness(val0[i_en],val_ref[i_en]);
   msg_log += sprintf("    -> %s compared rrad with no-rrad with goodness value: %.3e\n", ff, gn);
   
   if ( qualifier_exists("pl")){
      hplot(lo,hi,val0);
      ohplot(lo,hi,val_ref);
   }

   return gn;
}

define runtest(ffs){

   variable status = EXIT_SUCCESS;

   variable ff  = "relline_lp";
   
   if ( single_comparison_rrad(ff) < goodness_lim){
      msg_log += sprintf(" *** error: rrad line profile does not differ from standard profile\n");
      status = EXIT_FAILURE;
   }
   
   return status;
}
