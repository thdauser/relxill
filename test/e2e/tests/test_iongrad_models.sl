require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";

define check_single_beta(ff_ion){ %{{{

   
   variable lo,hi;
   variable val0, val1;
   (lo,hi) = log_grid(0.5,500,4000);

   fit_fun_default(ff_ion);

   set_par("*.beta",0);
   val0 = eval_fun_keV (lo,hi);
   
   
   set_par("*.beta",0.1);
   val1 = eval_fun_keV (lo,hi);
   
   

   variable gn = goodness(val0,val1);
   msg_log += sprintf("    -> %s: difference if beta changes? %.3e\n",
	    ff_ion, gn);
   return gn;
}
%}}}



define check_reflfrac_beta(ff_ion){ %{{{
   
   variable lo,hi;
   variable val0, val1;
   (lo,hi) = log_grid(0.5,500,4000);

   fit_fun_default(ff_ion);

   msg_log += sprintf(" checking the how refl_frac depends on BETA for %s\n",ff_ion);
   
   %% just use the normal relxilllp model for the test here
   set_par("*.ion_grad_type",0);
   set_par("*.boost",1,1);
   
   variable bet = [0:0.5:#3];
   variable ii, n = length(bet);

   _for ii(0,n-1){
      set_par("*.beta",bet[ii]);
      val0 = eval_fun_keV (lo,hi);
   }

   variable gn=0;
   
   return gn;
}
%}}}



define check_single_iongrad_model(ff_ion,ff_ref){ %{{{

   variable lo,hi;
   variable val0, val1, val2;
   (lo,hi) = log_grid(0.5,500,4000);

   putenv("RELXILL_NUM_RZONES=50");

   fit_fun_default(ff_ion);   
   set_par("*.xi_index",0);
   set_par("*.ion_grad_type",1);  %% PL ion gradient
   set_par("*.logN",15.0);
   val0 = eval_fun_keV (lo,hi);


   fit_fun_default(ff_ref);
   set_par("*.logN",15.0);
   val1 = eval_fun_keV (lo,hi);


   fit_fun_default(ff_ion);   
   set_par("*.ion_grad_type",2) ;  %% alpha disk
   val2 = eval_fun_keV (lo,hi);
   
   
   simple_plot(lo,hi,val0,val1;;__qualifiers());

   
   variable gn = goodness(val0,val1);
   msg_log += sprintf("    -> %s: goodness for comparison with relxilllp: %.3e\n",
	    ff_ion, gn);

   variable gn_change = goodness(val0,val2);
   msg_log += sprintf("    -> %s: difference between ion to non-ion grad? %.3e\n",
	    ff_ion, gn_change);

   putenv("RELXILL_NUM_RZONES");
   
   return gn;
}
%}}}


define runtest(ffs){ 

#ifndef STABLE

   counter++;
   msg_log += sprintf(" testing IONIZATION GRADIENT MODELS\n");


   %% currently the only model with an ionization gradient
   variable ff_ion = ["relxilllpionCp"];
   variable ff_ref = ["relxilllpCp"];
   
   variable ii,n = length(ff_ion);
   
   _for ii(0,n-1,1){
      if (check_single_iongrad_model(ff_ion[ii],ff_ref[ii];nopl) > goodness_lim*0.1){
	 msg_log += sprintf(" *** error: there seems to be a general problem with the IONGRAD MODEL %s \n",ff_ion[ii]);
	 return EXIT_FAILURE;
      }
      if (check_single_beta(ff_ion[ii];nopl) < goodness_lim*0.1){
	 msg_log += sprintf(" *** error: changing BETA paramter in IONGRAD MODEL seems to have no effect %s \n",ff_ion[ii]);
	 return EXIT_FAILURE;
      }
      if (check_reflfrac_beta(ff_ion[ii];nopl) > goodness_lim*0.1){
	 msg_log += sprintf(" *** error: reflection fraction in IONGRAD MODEL not working correctly %s \n",ff_ion[ii]);
	 return EXIT_FAILURE;
      }
   }   

#endif   
   
   return EXIT_SUCCESS;
   
}
%}}}

