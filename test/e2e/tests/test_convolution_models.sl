% -*- mode: slang; mode: fold -*-
require("load_test_setup.sl","Global");

variable msg_log = "";


define check_conv_mod_single(ff,ff_conv){ %{{{
   
   variable ener = 6.4;
   
   
   fit_fun_default(ff);
   set_par("*.lineE",ener);

   variable lo, hi;
   (lo,hi) = log_grid(0.2,10,600);
   
   variable val0 = eval_fun_keV(lo,hi);

   fit_fun(ff_conv+"(1,egauss)");
   set_par("*.center",ener);
   variable val1 = eval_fun_keV(lo,hi);
   
   variable sval1 = sum(val1);
   val1 *= sum(val0)/sval1;
   
   variable i_en = where(2.0<lo<6.0);
   simple_plot(lo[i_en],hi[i_en],val0[i_en],val1[i_en];;__qualifiers());
   
   variable gn = goodness(val0[i_en],val1[i_en]);
   msg_log += sprintf("    -> %s goodness value: %.3e\n",
	    ff, gn);
   return gn;
}
%}}}



define check_conv_norm_single(ff_conv){ %{{{
      
   variable ff_line = "egauss";
   fit_fun(ff_line);

   set_par("*.center",6.4);
   variable lo, hi;
   (lo,hi) = log_grid(0.2,10,600);
   
   variable val0 = eval_fun_keV(lo,hi);
      
   fit_fun(sprintf("%s(1,%s)",ff_conv,ff_line));
   variable val1 = eval_fun_keV(lo,hi);
         
   variable gn = sum(val0)/sum(val1); 
   msg_log += sprintf("    -> %s normalization value: %e\n",
	    ff_conv, gn);
   return gn;
}
%}}}


define runtest(ffs){ 
   
   counter++;
   msg_log += sprintf(" testing CONVOLUTION \n");
   
   variable ff = ["relline","relline_lp"];
   variable ff_conv = ["relconv","relconv_lp"];
   
   variable ii,n = length(ff);
   
   _for ii(0,n-1,1){
      if (not  (check_conv_mod_single(ff[ii],ff_conv[ii];nopl) < goodness_lim*8) ){
	msg_log += sprintf(" *** error: there seems to be a problem with the CONVOLUTION MODEL %s \n",ff_conv[ii]);
	 return EXIT_FAILURE;
      }
      
      if (not  ( abs(check_conv_norm_single(ff_conv[ii];nopl)-1) < 0.001) ){
	 msg_log += sprintf(" *** error: there seems to be a problem with the NORMALIZATION of the CONVOLUTION MODEL %s\n",ff_conv[ii]);
	 return EXIT_FAILURE;
      }

   }


   return EXIT_SUCCESS;
}

