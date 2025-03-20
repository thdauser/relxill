require("isisscripts");
require("tikz");

load_xspec_local_models("./build");

require("subs_model_plots.sl");


variable neval = 20;

variable ff = "relxilllpCp";

%%%
% overloads the basic function
variable global_refl_frac = 0;
define init_fitfun(_fitfun){
%%%
   fit_fun_default(_fitfun);
   set_par("*switch_reflfrac_boost",1);
   set_par("*.refl_frac",global_refl_frac,0,-1,10);
}

variable params = {};
% {pname, pmin, pmax, log}
list_append( params, {"a", 0,0.998,0});
list_append( params, {"h", 3,30,1} );
list_append( params, {"Incl", 20,70,0});
list_append( params, {"beta", 0,0.4,0});
list_append( params, {"Afe", 0.5,5,0});
% list_append( params, {"logxi", 1,4,0});
list_append( params, {"logN", 15,19,0});
%list_append( params, {"gamma", 1.8,2.6,0});



define param_multi_plot(_ff, _params){
   
   variable d_arr = {};
   variable p;
   foreach p(_params)
   {
      list_append( d_arr, get_model_evaluations(ff, p[0], p[1], p[2] ;
      log=p[3], neval=qualifier("neval",20)));
   }
   

   variable ii, nplots = length(d_arr);   
   variable nsub = int(ceil(nplots/2.0));

   if (qualifier_exists("single"))
   {
      variable pl_list = {};
      _for ii(0, nplots-1){
	 list_append(pl_list,plot_model_evaluations(d_arr[ii];noxtics=((ii==nplots-1)?0:1)));
      }
      return tikz_new_vbox_compound(__push_list(pl_list); interleave);
   }
   else
   {
      
      variable pl_list1 = {};
      _for ii(0, nsub-1){
	 list_append(pl_list1, plot_model_evaluations(d_arr[ii]; noxtics=((ii==nsub-1)?0:1))   );
      }
      
      variable pl_list2 = {};
      _for ii(nsub, nplots-1){
	 list_append(pl_list2,plot_model_evaluations(d_arr[ii];noxtics=((ii==nplots-1)?0:1)));
      }
      
      return tikz_new_hbox_compound(
      tikz_new_vbox_compound(__push_list(pl_list1); interleave),
      tikz_new_vbox_compound(__push_list(pl_list2); interleave),
      3; interleave );
   }
}

%putenv("RELXILL_RENORMALIZE=1");
%param_multi_plot(ff,neval).render("fig_relxill_renormalized.pdf");

putenv("RELXILL_RENORMALIZE=0");
variable global_refl_frac = -1;
variable rel_refl = param_multi_plot(ff,params; neval=10, single);
variable global_refl_frac = 0;
variable rel_prim = param_multi_plot(ff,params; neval=10, single);

tikz_new_hbox_compound(rel_refl, rel_prim, 3; interleave).render("fig_relxill_default_all.pdf");

variable global_refl_frac = 1;
param_multi_plot(ff,params; neval=10).render("fig_relxill_default.pdf");
