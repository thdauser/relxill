#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-


variable modlib = "build/librelxill.so";

if (stat_file(modlib) == NULL){
   message("\n **** error : local relxill model not found ; exiting ... **** \n ");
   exit;
}

require("xspec");
load_xspec_local_models(modlib);
require("isisscripts");

variable dir = "refdat/";

_traceback=1;

variable lo0, hi0;
(lo0,hi0) = log_grid(0.5,1000,3000);

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
   assoc["kTbb" ]   = "kTbb";
   assoc["frac_pl_bb" ]   = "Frac";
   assoc["A_CO" ]   = "A_CO";

   
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

putenv("DEBUG_RELXILL=1");

define check_xilltab_implementation_single(ff,tabname){
   variable tablepath =  getenv("RELXILL_TABLE_PATH")+ "/";
   
   variable ff_tab = "tab";
   
   add_atable_model(tablepath+tabname,"tab");
   
   fit_fun(ff);
   set_par("*.refl_frac",-1.0);

   if (string_matches(ff,"NS")!=NULL){
      set_par("*.kTbb",2.1);      
   } else {
      set_par("*.gamma",2.1);
   }
 
   set_par("*.Incl",40.0);
   
   variable val1     =  eval_fun_keV(lo0,hi0);
   list_par;
   
   variable pars = get_params();
   
   fit_fun(ff_tab);
   set_params_xillver(pars);
   variable valr =  eval_fun_keV(lo0,hi0);
   valr *=  sum(val1) / sum(valr);
   
   fits_write_binary_table(sprintf("%srefdat_%s.fits",dir,ff),
			   "REFDATA",
			  struct{lo=lo0,hi=hi0,val=valr}
			  );
   
   return ;

}


variable ff = __argv[1];
variable ff_tab = __argv[2];

check_xilltab_implementation_single(ff,ff_tab)