require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";


define check_neg_nthcomp_mod_single(ff,ff_dens){ %{{{
   
   
   variable lo,hi;
   variable val0, val1;
   (lo,hi) = log_grid(0.5,500,4000);
   
   fit_fun_default(ff_dens);
   set_par("*.kTe",100.0);
   val0 = eval_fun_keV (lo,hi);
      
   fit_fun_default(ff);
   set_par("*.Ecut",300.0);
   val1 = eval_fun_keV (lo,hi);
   
   simple_plot(lo,hi,val0,val1;;__qualifiers());

   variable gn = goodness(val0,val1);
   msg_log += sprintf("    -> %s vs %s  [ %.2e]\n",
	    ff_dens, ff,gn);
   return gn;
}
%}}}


define runtest(ffs){ 

   
   counter++;
   msg_log += sprintf(" testing the NTHCOMP and. BKN PL MODELS differ\n");
   
   variable ff = ["xillver","relxill","relxilllp"];
   variable ff_dens = ["xillverCp","relxillCp","relxilllpCp"];
   
   variable ii,n = length(ff);
   
   _for ii(0,n-1,1){
      if (check_neg_nthcomp_mod_single(ff[ii],ff_dens[ii];nopl) < goodness_lim*0.1){
	 msg_log += sprintf(" *** error: NTHCOMP and PKN PL do not differ \n => there seems to be a problem with the NTHCOMP MODEL %s \n",ff_dens[ii]);
	 return EXIT_FAILURE;
      }
   }

   

   return EXIT_SUCCESS;
   
}
