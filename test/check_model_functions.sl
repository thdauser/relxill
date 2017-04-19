#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-

_traceback=1;

%% we can load any model
variable ff = "*";
if (length(__argv)>1){
   ff = __argv[1];
}
load_xspec_local_models("build/librelxill.so");
require("isisscripts");

__set_hard_limits("relxilllp","h",-100,1000);


variable ALL_FF = ["relline","relline_lp","relxill","relxilllp","xillver","relxillD","xillverD"];
variable DATA_DIR = "refdata/";
variable goodness_lim = 1e-4;
variable sum_name = "plot_check_model_functions_sum.pdf";
variable counter = 0;

variable lo0, hi0;
(lo0,hi0) = log_grid(0.5,1000,3000);

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

   set_par("*.z",0.0);

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
   vmessage("    -> %s goodness value: %.3e",
	    ff, gn);
   return gn;
}
%}}}


define check_dens_mod_single(ff,ff_dens){ %{{{
   
   variable lxi = 2.0;
   
   variable lo,hi;
   variable val0, val1;
   (lo,hi) = log_grid(0.5,500,4000);
   
   fit_fun(ff_dens);
   set_par("*.refl_frac",1,0,-10,100);
   set_par("*.logN",15);
   set_par("*.logxi",lxi);
   val0 = eval_fun_keV (lo,hi);
   
   
   fit_fun(ff);
   set_par("*.refl_frac",1,0,-10,100);
   set_par("*.logxi",lxi);
   set_par("*.Ecut",300.0);
   val1 = eval_fun_keV (lo,hi);
   
   simple_plot(lo,hi,val0,val1;;__qualifiers());
   
   variable gn = goodness(val0,val1);
   vmessage("    -> %s goodness value: %.3e",
	    ff_dens, gn);
   return gn;
}
%}}}


define check_z(){ %{{{

   counter++;
   vmessage("\n### %i #### testing REDSHIFT parameter: ###",counter);
   
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
   vmessage("\n### %i ### testing SIMPLE EVALUATION ### ", counter);

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
   vmessage("\n### %i ### testing LINEE parameter: ###",counter);
   
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
   vmessage("\n### %i ### testing CONVOLUTION: ###",counter);
   
   variable ff = ["relline"];
   variable ff_conv = ["relconv"];
   
   variable ii,n = length(ff);
   
   _for ii(0,n-1,1){
      if (check_conv_mod_single(ff[ii],ff_conv[ii];nopl) > goodness_lim*5){
	vmessage(" *** error: there seems to be a problem with the CONVOLUTION MODEL %s ",ff_conv[ii]);
	 return EXIT_FAILURE;
      }
   }


   return EXIT_SUCCESS;
}
%}}}

define check_dens_mod(){ %{{{
   
   counter++;
   vmessage("\n### %i ### testing HIGH DENSITY MODELS: ###",counter);
   
   variable ff = ["xillver","relxill"];
   variable ff_dens = ["xillverD","relxillD"];
   
   variable ii,n = length(ff);
   
   _for ii(0,n-1,1){
      if (check_dens_mod_single(ff[ii],ff_dens[ii];nopl) > goodness_lim*0.1){
	vmessage(" *** error: there seems to be a problem with the HIGH DENSITY MODEL %s ",ff_dens[ii]);
	 return EXIT_FAILURE;
      }
   }

   

   return EXIT_SUCCESS;
   
}
%}}}

define test_caching_spec(v,v0,inp,inp0){ %{{{
   variable ii, n = length(v);
   variable status = EXIT_SUCCESS;
   
   variable lim = 1e-8;
   variable tmp;
   
   % #1# for v each spectrum should be different
   _for ii(0,n-2,1){
      tmp = sum(  (v[ii] - v[ii+1] )^2 ) ;
      if (tmp < lim ){
	 vmessage("     caching spectra not correct (did not change while it should), %.3e < %.3e for values %.2e and %.2e ",
		  tmp, lim, inp[ii], inp[ii+1] );
	 return EXIT_FAILURE;
      }
      
   }
   
   % #2# for v0 each spectrum should be the same
   _for ii(0,n-2,1){
      tmp = sum(  (v0[ii] - v0[ii+1] )^2 );
      if (tmp > lim){
	 vmessage(" *** caching spectra not correct (did change while it should NOT), %.3e > %.3e", tmp, lim);
	 return EXIT_FAILURE;
      }
      
   }

   % #3# by design the last elements should be the same
   tmp = sum(  (v[-1] - v0[-1] )^2 );
   if (tmp > lim){
      vmessage(" *** caching spectra not correct (spectra are different while they should NOT), %.3e > %.3e", tmp, lim);
      return EXIT_FAILURE;
   }
   
   
   return EXIT_SUCCESS;
}
%}}}

define check_caching_single(ff,par){ %{{{

   fit_fun(ff);   
   variable param0 = get_params("*."+par);
   if (length(param0)!=1){
      vmessage(" *** error *** problem with model %s and parameter %s when testing caching",ff,par);
      return EXIT_FAILURE;
   }
   variable p = param0[0];
   % clone parameter

   % ### 1 ### change parameter between min and max
   variable N = 3;
   variable v;
   
   variable vals;

   %% make sure we get useful values of p.value==0
   if (p.value < 0.1 && p.value >= 0){
      p.value = 0.1;
   }
   
   if (abs(p.value-p.min) < abs(p.max-p.value)){
      vals = [p.value:p.value*1.1 :#N];
   } else {
      if (p.value>=0){
	 vals = [p.value*0.9:p.value:#N];
      } else {
	 vals = [p.value:p.value*1.1 :#N];	 
      }
   }

   variable vals0 = vals*0 + vals[-1];  %% DO NOT CHANGE THIS (will be used later)
   
   %% make sure values of the table are loaded before such that we do
   %% not count the loading of the table as well here
   foreach v(vals){
      set_par(p.name,v);
      () = eval_fun_keV(lo0,hi0);
   }
   
   variable val_ncache = Array_Type[N];				      
   variable val0_cache = Array_Type[N];				      

   %% special cases %%
   if (string_match(p.name,`.*\.Rbr`)){
      set_par("*.Index1",3.5);
   }
   %% %%%%%%%%%%%%% %%
   
   variable ii;
   tic;
   _for ii(0,N-1){
      set_par(p.name,vals[ii]);
      val_ncache[ii] = eval_fun_keV(lo0,hi0);
   }
   variable dt = toc;

   tic;
   _for ii(0,N-1){
      set_par(p.name,vals0[ii]);
      val0_cache[ii] = eval_fun_keV(lo0,hi0);
   }
   variable dt2 = toc;

   if (test_caching_spec(val_ncache,val0_cache,vals,vals0) != EXIT_SUCCESS){
      vmessage("     => spectra did not behave as expected when testing caching (%s -  %s)",ff,par);
      return EXIT_FAILURE;
   } else {
      if (dt2/N*1e3 < 25.0){
	 vmessage("   --> %.1f ms vs. %.1f ms  --  (%s -  %s)",dt/N*1e3,dt2/N*1e3,ff,par );	 
      } else {	 
	 vmessage("   --> *** FAILED ***  [ %.1f ms vs. %.1f ms ] --  (%s -  %s) ",dt/N*1e3,dt2/N*1e3,ff,par );	       
	 return EXIT_FAILURE;
      }
   }
   
   
   return EXIT_SUCCESS;
}
%}}}

define check_caching(){ %{{{

   counter++;
   vmessage("\n### %i ### testing CACHING for following models ###",counter);
   
   variable ff_arr = Assoc_Type[Array_Type];
   
   variable std_rel_param = ["a","Incl","Rin","Rout"];
   variable std_xill_param = ["logxi","Afe","z"];
   
   ff_arr["relxill"]   = [std_rel_param, "Rbr" , "Index1","Index2",std_xill_param, "Ecut"];
   ff_arr["relxilllp"] = [std_rel_param, "h","refl_frac", std_xill_param, "Ecut" ];
   ff_arr["relxillD"]   = [std_rel_param, "Rbr", "Index1","Index2", std_xill_param, "logN" ];
   
   variable ff, params;
   variable ii, n;
   foreach ff(assoc_get_keys(ff_arr)){
      params = ff_arr[ff];
      n = length(params);
      _for ii(0,n-1){
	 if (check_caching_single(ff,params[ii];nopl) == EXIT_FAILURE){
	    vmessage(" *** error: there seems to be a problem with the CACHING for model %s and parameter %s",
		     ff,ff_arr[ff][ii]);
	    return EXIT_FAILURE;
	 }
      }      
   }
   return EXIT_SUCCESS;   
}
%}}}

define do_mc_testing(){ %{{{

   
   counter++;
   vmessage("\n### %i ### random MC Parameter tests ### ", counter);

   variable n_mc = 100;
   variable stat;

   variable lo,hi;
   (lo,hi) = log_grid(0.1,500,500);
   variable fakeval = lo^-(1.9);
   variable fake_dat = struct{
         bin_lo=lo,
         bin_hi=hi,
         value=fakeval,
         err= 0.1*fakeval
   };
   
   variable iDat = define_counts(_A(fake_dat));

   variable save_fit_verbose = Fit_Verbose;
   Fit_Verbose=0;
   
   
   variable ffs = ALL_FF;
   variable ff;
   variable val, dt;
   foreach ff(ffs){

      fit_fun(ff);
      freeze("*.norm");
      tic;
      stat = fit_search(n_mc, &eval_counts;  serial, dir="/tmp/fit_search");
      dt = toc / n_mc * 1e3;
      if ( stat == NULL ){
	 vmessage(" *** error: MC Parameter test for %s failed!",ff);
	 return EXIT_FAILURE;
      }
      vmessage("   -> successfully tested %s with %i evaluations   \t [ <time> ~ %.0fmsec ]",ff,n_mc,dt);

   }
   
   Fit_Verbose = save_fit_verbose;
   delete_data(iDat);
   
   return EXIT_SUCCESS;
}
%}}}


if (eval_test() != EXIT_SUCCESS) exit;

if (check_z() != EXIT_SUCCESS) exit;
if (check_linee() != EXIT_SUCCESS) exit;
if (check_conv_mod() != EXIT_SUCCESS) exit;
if (check_dens_mod() != EXIT_SUCCESS) exit;
if (check_caching() != EXIT_SUCCESS) exit;

if (do_mc_testing() != EXIT_SUCCESS) exit;
