require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_model.sl");

load_relxill_dev_model();

_traceback=1;

variable ff = "relxillBBret";
fit_fun_default(ff);

variable params = Assoc_Type[Array_Type];

params["a"] = [0.85,0.851,0.9,0.901];


%%%% FUNCTIONS %%%%% 


%% overwrite the general routine %%%
define plot_single_param(pars,key){
   
   variable param_vals = pars[key];
   variable ii, n = length(param_vals);
   
   variable col = [CB_COLOR_SCHEME_NB,CB_COLOR_SCHEME_NB];
   
   variable pl = get_std_plot();
   
   variable lab = String_Type[n];
   variable keyTexSafe = strreplace(key,"_","-");
   
   fit_fun_default(ff);
   variable egrid = get_ener_grid();
   _for ii(0,n-1){
           
      set_par("*.refl_frac",-1.0,0,-10,10);
      egrid = calc_single_evaluation(pl,egrid,key,param_vals[ii]);
      pl.plot(0.5*(egrid.lo+egrid.hi), egrid.val; color=col[ii],line=1);

      set_par("*.refl_frac",0.0,0,-10,10);
      egrid = calc_single_evaluation(pl,egrid,key,param_vals[ii]);
      pl.plot(0.5*(egrid.lo+egrid.hi), egrid.val; color=col[ii],line=1);

      set_par("*.refl_frac",1.0,0,-10,10);
      egrid = calc_single_evaluation(pl,egrid,key,param_vals[ii]);
      pl.plot(0.5*(egrid.lo+egrid.hi), egrid.val; color=col[ii],line=0);

      vmessage("a = %.3e : Integ Flux = %e ", param_vals[ii], sum(egrid.val) );
      
      lab[ii] = sprintf("%s = %.3e", keyTexSafe, param_vals[ii]);
   }
     
   add_plot_info(pl,lab,keyTexSafe,col);
   
   return pl;
}



variable pl = plotAllParameterCombinations(params);
variable combpl = xfig_multiplot(pl);
combpl.render(dir_plots+"plotDifferentSpins.pdf");
