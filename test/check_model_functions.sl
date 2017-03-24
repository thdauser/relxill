% #!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-

%% we can load any model
variable ff = "*";
if (length(__argv)>1){
   ff = __argv[1];
}
load_xspec_local_models("build/librelxill.so");
require("isisscripts");

__set_hard_limits("relxilllp","h",-100,1000);


variable ALL_FF = ["relline","relline_lp","relxill","relxilllp","xillver"];
variable DATA_DIR = "refdata/";
variable goodness_lim = 1e-4;
variable sum_name = "plot_check_model_functions_sum.pdf";
variable counter = 0;


%%% generally useful functions %%% 

define get_single_name(){ %{{{
   counter++;
   return sprintf("check_model_functions_%04i.pdf");
}
%}}}
define goodness(d1,d2){ %{{{
   variable ind = where(d2!=0);
   return  sqrt(sum((d1[ind]/d2[ind]-1)^2))/length(d1[ind]);   
}
%}}}
define simple_plot(lo,hi,val0,val1){ %{{{
   if (qualifier_exists("nopl")){
      return;
   }
   
   open_plot("/xw",2,1);
   hplot(lo,hi,val0);
   ohplot(lo,hi,val1);
   
   variable ind = where(val1!=0);
   hplot(lo[ind],hi[ind],val0[ind]/val1[ind]);
}
%}}}


variable EXIT_SUCCESS = 0;
variable EXIT_FAILURE = -1;

variable counter = 0;

define check_z_param(ff){ %{{{
   fit_fun(ff);
   
   variable zp = 0.1;
   
   variable lo, hi;
   (lo,hi) = log_grid(0.2,900,3000);
   
   set_par("*.z",0);
   variable val = eval_fun_keV(lo,hi);
      
   set_par("*.z",zp);
   variable val_z = eval_fun_keV(lo,hi);

   % now change grid
   variable zlo = lo * (1+zp) ;
   variable zhi = hi * (1+zp) ;
   
   variable val_z_corr = rebin(lo,hi,zlo,zhi,val_z);


   
   
   variable i_en;
   
   if (string_match(ff,"\.*line") == 0){
      i_en = where(1.0 < lo < 120.0);
   } else {
      i_en = where(1.0 < lo < 6.0);
   }
	
   lo = lo[i_en];
   hi = hi[i_en];
   zlo = zlo[i_en];
   zhi = zhi[i_en];
   val = val[i_en];
   val_z = val_z[i_en];
   val_z_corr = val_z_corr[i_en];
   
   variable ind = where(val_z_corr!=0);

   variable fac = (0.5*(lo+hi))/(hi-lo);
   
   if (not (qualifier_exists("nopl"))) {
      xlog;
      ylog;
      open_plot("/xw",2,1);
      hplot(lo,hi,val*fac);
      ohplot(lo,hi,val_z_corr*fac);
%      ohplot(zlo,zhi,val_z/(zhi-zlo)*(0.5*(zlo+zhi))^2);
      ohplot(lo,hi,val_z*fac);
      
      ylin;
      yrange(0.95,1.05);
      hplot(lo[ind],hi[ind],val[ind]/val_z[ind]);
      ohplot(lo[ind],hi[ind],val[ind]/val_z_corr[ind]);
   }
   
   variable gn = goodness(val,val_z_corr);
   vmessage("    -> %s goodness value: %.3e",
	    ff, gn);
   return gn;
}
%}}}


define check_line_ener(ff){ %{{{
   fit_fun(ff);

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
   
   variable i_en = where(2.0<lo<6.0);
   variable gn = goodness(val0[i_en],val_reb[i_en]);
   vmessage("    -> %s goodness value: %.3e",
	    ff, gn);
   return gn;
}
%}}}


define check_conv_mod_single(ff,ff_conv){ %{{{
   
   variable ener = 6.4;
   
   
   fit_fun(ff);
   set_par("*.lineE",ener);
   list_par;
   variable lo, hi;
   (lo,hi) = log_grid(0.2,10,1200);
   
   variable val0 = eval_fun_keV(lo,hi);
      
   fit_fun(ff_conv+"(1,egauss)");
   set_par("*.center",ener);
   list_par;
   variable val1 = eval_fun_keV(lo,hi);
   
   variable sval1 = sum(val1);
   val1 *= sum(val0)/sval1;
   
   variable i_en = where(2.0<lo<6.5);
   simple_plot(lo[i_en],hi[i_en],val0[i_en],val1[i_en];;__qualifiers());
   
   variable gn = goodness(val0[i_en],val1[i_en]);
   vmessage("    -> %s goodness value: %.3e",
	    ff, gn);
   return gn;
}
%}}}


define check_z(){ %{{{

   counter++;
   vmessage("\n ### %i #### testing REDSHIFT parameter: ###",counter);
   
   variable ffs = ALL_FF;
   variable ff;
   
   foreach ff(ffs){            
      if (check_z_param(ff;nopl) > goodness_lim){
	vmessage(" *** error: there seems to be a problem with the redshift in %s ",ff);
	 return EXIT_FAILURE;
      }
   }

   
   message(" ");

   return EXIT_SUCCESS;
}
%}}}

define eval_test(){ %{{{

   counter++;
   vmessage("\n ### %i ### testing SIMPLE EVALUATION ### ", counter);

   variable ffs = ALL_FF;
   variable ff;
   variable val;
   foreach ff(ffs){
      fit_fun(ff);
      val = eval_fun_keV(1,2);
      if (not ( val > 0 )){
	 vmessage(" *** error: simple test for %s failed!",ff);
	 return EXIT_FAILURE;
      }
   }
   
   return EXIT_SUCCESS;
}
%}}}

define check_linee(){ %{{{
   
   counter++;
   vmessage("### %i ### testing LINEE parameter: ###",counter);
   
   variable ffs = ["relline","relline_lp"];
   variable ff;
   
   foreach ff(ffs){            
      if (check_line_ener(ff;nopl) > goodness_lim*3){
	vmessage(" *** error: there seems to be a problem with the LineE in %s ",ff);
	 return EXIT_FAILURE;
      }
   }

   
   message(" ");

   return EXIT_SUCCESS;
}
%}}}


define check_conv_mod(){ %{{{
   
   counter++;
   vmessage("### %i ### testing CONVOLUTION: ###",counter);
   
   variable ff = ["relline"];
   variable ff_conv = ["relconv"];
   
   variable ii,n = length(ff);
   
   _for ii(0,n-1,1){
      if (check_conv_mod_single(ff[ii],ff_conv[ii]) > goodness_lim*5){
	vmessage(" *** error: there seems to be a problem with the CONVOLUTION MODEL %s ",ff_conv[ii]);
%	 return EXIT_FAILURE;
      }
   }

   
   message(" ");

   return EXIT_SUCCESS;
}
%}}}



if (eval_test() != EXIT_SUCCESS) exit;
if (check_z() != EXIT_SUCCESS) exit;
if (check_linee() != EXIT_SUCCESS) exit;
if (check_conv_mod() != EXIT_SUCCESS) exit;

