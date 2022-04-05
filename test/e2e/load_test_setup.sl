#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-

require("subs_fitfunctions.sl");
variable modlib = "../build/librelxill.so";
load_relxill_model_devel(modlib);

require("isisscripts");

variable counter = 0;
_traceback=1;


variable goodness_lim = 1e-4;

variable EXIT_SUCCESS = 1;
variable EXIT_FAILURE = 0;

variable msg_log = NULL;
 
variable sum_name = "plot_check_model_functions_sum.pdf";

variable lo0, hi0;
(lo0,hi0) = log_grid(0.1,1000,3000);


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

   xlog;ylog;
   open_plot("/xw",2,1);
   hplot(lo,hi,val0);
   ohplot(lo,hi,val1);
   
   variable ind = where(val1!=0);
   ylin;
   hplot(lo[ind],hi[ind],val0[ind]/val1[ind]);
}
%}}}
define grav_redshift_prim(a,h){ %{{{
   return 1.0 / sqrt( 1.0 - 2*h/(h^2 + a^2) ) - 1.0;
}
%}}}
define set_params_xillver(pars){ %{{{

   if (pars==NULL) return;
   variable ii, n = length(pars);

   variable assoc = Assoc_Type[String_Type];
   assoc["gamma"]   = "Gamma";
   assoc["Afe"]     = "A_Fe";
   assoc["logxi"]   = "logXi";
   assoc["norm"]    = "norm";
   assoc["Ecut"]    = "Ecut";
   assoc["Incl"]    = "Incl";
   assoc["z"   ]    = "redshift";
   assoc["logN"]    = "Dens";
   assoc["kTe" ]    = "kTe";
   assoc["frac_pl_bb"]     = "Frac";
   assoc["kTbb"]     = "kTBB";

   
   variable par_array = ["gamma","Afe"];
   
   _for ii(0,n-1){
      %% get the name of the parameter string and then set the
      %% respective value
      variable pname = string_matches(pars[ii].name,"\.\(.*\)"R)[-1];
      if (pname == "refl_frac") continue; 
      
      set_par("*"+assoc[pname],pars[ii].value);
   }
   
   
}
%}}}
define isLpModel(ff){ %{{{
   if (string_matches(ff,"relxilllp")!=NULL){
      return 1;
   } else {
      return 0;
   }
}
%}}}


define is_lamppost_model(ff){ %{{{
   if (string_matches(ff,"lp")!=NULL){
      return 1;
   } else {
      return 0;
   }
}
%}}}


define is_conv_model(ff){ %{{{
   if (string_matches(ff,"conv")!=NULL){
      return 1;
   } else {
      return 0;
   }
}
%}}}


define is_xillver_model(ff){ %{{{
   if (string_matches(ff,"xillver")!=NULL){
      return 1;
   } else {
      return 0;
   }
}
%}}}

define is_relxill_model(ff){ %{{{
   if (string_matches(ff,"relxill")!=NULL){
      return 1;
   } else {
      return 0;
   }
}
%}}}

define is_line_model(ff){ %{{{
   if (string_matches(ff,"line")!=NULL){
      return 1;
   } else {
      return 0;
   }
}
%}}}

define model_exists(ff,ffs){ %{{{
   if (length(where(ff == ffs )) ==0 ){
      return 0;
   } else {
      return 1;
   }
}

%}}}
