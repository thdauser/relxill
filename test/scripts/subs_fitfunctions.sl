

define load_relxill_model_devel(modlib){
   
   if (stat_file(modlib) == NULL){
      message("\n **** error : local relxill model not found ; exiting ... **** \n ");
      exit(1);
   }
   
   require("xspec");
   load_xspec_local_models(modlib);
   
   
   __set_hard_limits("relxilllp","h",-100,1000);

   fit_fun("relxill");
   () = eval_fun(1,2);
}


variable ALL_FF;  %% empty global variable

define get_implemented_fitfunctions(){
   variable ALL_FF = ["relline","relline_lp","relxill","relxilllp","xillver",
		  "relxillCp","relxilllpCp","xillverCp","relxilllpion","relxilllpionCp"];
   variable additional_FF =  ["xillverNS","relxillNS","xillverCO","relxillCO","relxillDCp","relxilllpDCp","xillverDCp"];
   
   if (qualifier("dev",0) == 1 ){
      ALL_FF = [ALL_FF, additional_FF];
   }
   if (qualifier_exists("only_dev") ){
      return additional_FF;
   }
   return ALL_FF;
}
   


define isFitfunDevModel(ff){ %{{{
   if (length( where( ff==get_implemented_fitfunctions(;dev=1) )) >0){
      return 1;
   } else {
      return 0;
   }
}
%}}}


define fit_fun_default(ff){ %{{{
   %% only works for single model components %%
   
   fit_fun(ff);
   
   variable p,pa = get_params();   
   variable df;
   foreach p(pa){
      df = eval(sprintf("%s_default(%i)",qualifier("ff0",ff),p.index-1));
      p.value=df.value;
      p.min=df.min;
      p.max=df.max;
      p.freeze=df.freeze;
   }
   set_params(pa);
}
%}}}


define randomize_params() { %{{{
   variable p = get_params();
   
   variable i;   
   for (i = 0; i < length(p); i++)   {
      if (p[i].freeze == 0 and p[i].tie == NULL) {
	 p[i].value = rand_flat (p[i].min, p[i].max);
      }
   }
   
   set_params(p);
}
%}}}




define selectModelByType(ff, strType){ %{{{
   variable ii, n = length(ff);
   variable selectedFf = String_Type[0];
   _for ii(0, n-1){
      if (string_match(ff[ii], strType) != 0){
	 selectedFf = [selectedFf, ff[ii]];
      }
   }
   return selectedFf;
}
%}}}
define getAllRelxillTypeModels(){ %{{{
   return selectModelByType(ALL_FF, "relxill");
}
%}}}
define getAllXillverTypeModels(){ %{{{
   return selectModelByType(ALL_FF, "xillver");   
}
%}}}
