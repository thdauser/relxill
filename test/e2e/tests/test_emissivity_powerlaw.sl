require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";

load_relxill_model_devel(modlib);


define test_linemodel_normalization(){ %{{{
   
   variable status = EXIT_SUCCESS;
   

   
   fit_fun("relline");
   
   set_par("*.Rbr",400);
   set_par("*.Rout",1000,1,1,1000);
   variable rin = [-1,-100];
   
   set_par("*.Rin",-1);
   variable val1 = eval_fun_keV(0.1,10);
   set_par("*.Rin",-50);
   variable val2 = eval_fun_keV(0.1,10);

   msg_log += sprintf(" - relline normalization is %.3e for Rin=R_ISCO and %.3e for Rin=100*R_ISCO\n", val1, val2);
   
   
   
   if (val1 >= val2*0.9){
      vmessage(" *** error: relline for Rin=R_ISCO does not have lower normalization as truncated disk");
      status = EXIT_FAILURE;
   }

   if (val1 >= 1.0){
      vmessage(" *** error: relline for Rin=R_ISCO must have lower normalization as 1 (value is %e)",val1);
      status = EXIT_FAILURE;
   }

   variable norm_expect = cos(30./180*PI)*0.5;
   if ( abs(val2-norm_expect) > 0.05 ){
      vmessage(" *** error: relline norm for highly truncated disk must be close to 0.5*cos(incl)=%.3f (value is %.4f)",norm_expect, val2);
      status = EXIT_FAILURE;
   }



   return status;
}
%}}}


define test_convmodel_normalization(){ %{{{
   
   variable status = EXIT_SUCCESS;
   
   variable ff_r = "relconv(1,powerlaw)";
   variable ff_p = "powerlaw";
   

   variable lo0,hi0;
   (lo0, hi0) = log_grid(0.1,500,1000);
   
   fit_fun(ff_r);
   
   set_par("*.PhoIndex",2);
   set_par("*.Rbr",400);
   set_par("*.Rout",1000,1,1,1000);
   variable rin = [-1,-100];
   
   set_par("*.Rin",-1);
   variable val1 = eval_fun_keV(lo0,hi0);
   fit_fun(ff_p);
   variable val1p = eval_fun_keV(lo0,hi0);
   

   
   
   
   fit_fun(ff_r);
   set_par("*.Rin",-100);
   variable val2 = eval_fun_keV(lo0, hi0);
   fit_fun(ff_p);
   variable val2p = eval_fun_keV(lo0,hi0);
   
   variable ind = where(0.1 < lo0 < 500.0);
   
   variable rat_risco = sum(val1[ind]) / sum(val1p[ind]);
   variable rat_trunc = sum(val2[ind]) / sum(val2p[ind]);
   
   
   msg_log += sprintf(" - relconv normalization is %.3e for Rin=R_ISCO and %.3e for Rin=100*R_ISCO\n", rat_risco, rat_trunc);
   
   
   
   if (rat_risco >= rat_trunc){
      vmessage(" *** error: relconv for Rin=R_ISCO does not have lower normalizationt as truncated disk");
      status = EXIT_FAILURE;
   }

   if (rat_risco >= 1.0){
      vmessage(" *** error: relconv for Rin=R_ISCO must have lower normalization as 1 (value is %e)",rat_risco);
      status = EXIT_FAILURE;
   }

   variable norm_expect = cos(30./180*PI)*0.5;
   if ( abs(rat_trunc-norm_expect) > 5e-2 ){
      vmessage(" *** error: relconv for highly truncated disk must be close to 0.5*cos(incl)=%.3f (value is %.4f)",norm_expect, rat_trunc);
      status = EXIT_FAILURE;
   }


   
   return status;
}
%}}}




define runtest(ffs){ %{{{
      
   variable status = EXIT_SUCCESS; 

   variable ff;
   putenv("RELLINE_PHYSICAL_NORM=1");

   variable local_status = EXIT_SUCCESS;

   variable status_line = test_linemodel_normalization();
   variable status_conv = test_convmodel_normalization();	 

   
   if (( status_line == EXIT_FAILURE ) || (status_conv == EXIT_FAILURE)) {
      status = EXIT_FAILURE;
      vmessage(" *** error: test failed ");
   }
   
   
   putenv("RELLINE_PHYSICAL_NORM");

   
   return status;
}
%}}}


