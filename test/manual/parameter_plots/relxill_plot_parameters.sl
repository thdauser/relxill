require("isisscripts");
require("tikz");

load_xspec_local_models("./build");

require("./subs_model_plots.sl");


variable neval = 20;

variable ff = "relxilllpCp_ext";

%%%
% overloads the basic function
variable global_refl_frac = 0;
define init_fitfun(_fitfun){
%%%
   __set_hard_limits (_fitfun, "Rin", -100, 200; parm=[-100, -1, 100]);
   fit_fun_default(_fitfun);
   set_par("*switch_reflfrac_boost",1);
   set_par("*.refl_frac",global_refl_frac,0,-1,10);
}

variable params = {};
% {pname, pmin, pmax, log}
list_append( params, {"a", 0,0.998,0});
list_append( params, {"h", 3,30,1} );
list_append( params, {"Incl", 20,70,0});
% list_append( params, {"beta", 0,0.4,0});
list_append( params, {"r_ring", 0.0,10.0,0});
list_append( params, {"Afe", 0.5,5,0});
% list_append( params, {"logxi", 1,4,0});
list_append( params, {"logN", 15,19,0});
list_append( params, {"gamma", 1.8,2.6,0});
list_append( params, {"Rin", 1.3,30,0});



%putenv("RELXILL_RENORMALIZE=1");
%param_multi_plot(ff,neval).render("fig_relxill_renormalized.pdf");

putenv("RELXILL_RENORMALIZE=0");
%variable global_refl_frac = -1;
%variable rel_refl = param_multi_plot(ff,params; neval=10, single);
%variable global_refl_frac = 0;
%variable rel_prim = param_multi_plot(ff,params; neval=10, single);

%tikz_new_hbox_compound(rel_refl, rel_prim, 3; interleave).render("fig_relxill_default_all.pdf");

variable global_refl_frac = 1;
variable plt = param_multi_plot(ff,params; neval=10);

print(plt);
plt.render("fig_relxill_default.pdf");
