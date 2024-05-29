require("isisscripts");
require("tikz");

load_xspec_local_models("./build");

require("subs_model_plots.sl");


variable neval = 20;

putenv("RELXILL_RENORMALIZE=1");
variable ff = "relxilllpCp";

variable d_arr = {};
list_append( d_arr, get_model_evaluations(ff, "a", 0,0.998; neval=neval));
list_append( d_arr, get_model_evaluations(ff, "h", 3,30; log, neval=neval));
list_append( d_arr, get_model_evaluations(ff, "Incl", 20,70; neval=neval));
list_append( d_arr, get_model_evaluations(ff, "refl_frac", 0.5,2; neval=neval));
list_append( d_arr, get_model_evaluations(ff, "Afe", 0.5,5; neval=neval));
list_append( d_arr, get_model_evaluations(ff, "logxi", 1,4; neval=neval));
list_append( d_arr, get_model_evaluations(ff, "logN", 15,19; neval=neval));
list_append( d_arr, get_model_evaluations(ff, "gamma", 1.8,2.6; neval=neval));


variable d;
variable ii, nplots = length(d_arr);

variable nsub = int(ceil(nplots/2.0));

variable pl_list1 = {};
_for ii(0, nsub-1){
   list_append(pl_list1, plot_model_evaluations(d_arr[ii]; noxtics=((ii==nsub-1)?0:1))   );
}

variable pl_list2 = {};
_for ii(nsub, nplots-1){
   list_append(pl_list2,plot_model_evaluations(d_arr[ii];noxtics=((ii==nplots-1)?0:1)));
}

tikz_new_hbox_compound(
tikz_new_vbox_compound(__push_list(pl_list1); interleave),
tikz_new_vbox_compound(__push_list(pl_list2); interleave),
3; interleave ).render("fig_relxill_normalized.pdf");
