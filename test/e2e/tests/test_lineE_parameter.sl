% -*- mode: slang; mode: fold -*-
require("load_test_setup.sl","Global");

variable msg_log = "";


define check_line_ener(ff){ %{{{
   fit_fun_default(ff);

   if (string_match(ff,"\.*line") == 0){
      return;
   }
      
   variable le0 = 6.4;
   variable le1 = 3.2;
   
   variable lo, hi;
   (lo,hi) = log_grid(0.2,10,600);
   
   set_par("*.lineE",le0);
   variable val0 = eval_fun_keV(lo,hi);
      
   set_par("*.lineE",le1);
   variable val1 = eval_fun_keV(lo,hi);

   % now change grid
   variable lo1 = lo * (le0/le1) ;
   variable hi1 = hi * (le0/le1) ;
   
   variable val_reb = rebin(lo,hi,lo1,hi1,val1);
   
   if (not qualifier_exists("nopl")){
      hplot(lo,hi,val0);
      ohplot(lo,hi,val_reb);
      sleep(5);      
   }
   
   variable i_en = where(2.0<lo<6.0);
   variable gn = goodness(val0[i_en],val_reb[i_en]);
   msg_log += sprintf("    -> %s goodness value: %.3e\n", ff, gn);
   return gn;
}
%}}}

define check_line_limb(ff){ %{{{
   fit_fun_default(ff);

   if (string_match(ff,"\.*line") == 0){
      return;
   }
      
   variable le0 = 6.4;
   variable le1 = 3.2;
   
   variable lo, hi;
   (lo,hi) = log_grid(0.2,10,600);
   
   
   variable val0 = eval_fun_keV(lo,hi);
      
   set_par("*.limb",1);
   variable val1 = eval_fun_keV(lo,hi);
   
   set_par("*.limb",2);
   variable val2 = eval_fun_keV(lo,hi);

   variable gn = goodness(val0,val1);
   variable gn2 = goodness(val0,val2);
   msg_log += sprintf("    -> %s goodness value for limb dark %.3e and bright %.3e\n",
	    ff, gn,gn2);
   return gn+gn2;
}
%}}}


define runtest(ffs){
   
   msg_log += sprintf(" testing LINEE parameter\n");

   variable ff;
   foreach ff(ffs){            
      
      if(is_line_model(ff)){
	 
	 if (check_line_ener(ff) > goodness_lim*5){
	    msg_log += sprintf(" *** error: there seems to be a problem with the LineE in %s \n",ff);
	    return EXIT_FAILURE;
	 }
	 if (check_line_limb(ff) < goodness_lim*5){
	    msg_log += sprintf(" *** error: there seems to be a problem with the Limb Brightening / Darkening  in %s \n",ff);
	    return EXIT_FAILURE;
	 }
      }
   }
      
   return EXIT_SUCCESS;
}
