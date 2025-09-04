require("isisscripts");

require("subs_jedsad.sl");

require("/home/thomas/git/relxill.git/test/manual/parameter_plots/subs_model_plots.sl");



variable neval_global = 5;

load_models();
setup_model();

variable param_list = ["Risco","Rj", "Mdot", "Omega", "Ms"];
variable is_log = [0,0,1,1,1];

define single_eval(){

   variable norm_energy = 10.0; %% keV
   
   fit_fun(global_fitfun);
   variable full = eval_ff();

   
   fit_fun("relxill_jedsad");
   set_par("*.refl_frac",-1);
   variable relxill_refl = eval_ff();

   fit_fun("relxill_jedsad");
   set_par("*.refl_frac",0);
   variable relxill_prim = eval_ff();

   variable norm_jed = relxill_prim.val[wherefirst(relxill_prim.lo > norm_energy)];
   

   fit_fun("jedsad");
   variable jed = eval_ff(;norm=norm_jed, norm_energy=norm_energy);

   %% reset fit function
   fit_fun(global_fitfun);

   return struct_combine(full,
   struct{jed=jed.val, relxill_refl=relxill_refl.val, relxill_prim=relxill_prim.val}
   );

   
   
}

define eval_model_for_param(pname){

   setup_model();
   
   variable pinfo = get_params("jedsad*.$pname"$)[0];      
   variable val_min = pinfo.min;
   variable val_max = pinfo.max;
   variable val = pinfo.value;

   variable neval = qualifier("neval", neval_global);
   variable dat = Struct_Type[neval];
   
   variable parvalues = get_parvalues(val_min, val_max, neval ;;__qualifiers() );
   variable full_params = Array_Type[neval];

   variable ii;
   _for ii(0, neval-1){
      vmessage(" Setting Paramter %s (%i) to %.3e", pname, pinfo.index, parvalues[ii]);
      set_par(pinfo.index, parvalues[ii]);
      dat[ii] = single_eval();
   }
      
   %% reset the parameter
   set_params(pinfo);

   return dat, parvalues;
}

define single_plot_data(pl, da,col,opa){

   pl.plot(da.emean, da.relxill_prim;
   opacity=opa, color=col, width=2);

   if (0 < sum(da.relxill_refl) < 1e10){
      pl.plot(da.emean, da.relxill_refl;
      opacity=opa, color=col, width=1);
   }

   pl.plot(da.emean, da.jed;
   opacity=opa, color=col, width=2, line=1);

}

define plot_data_struct(data, parvalues, parname){

   variable n = length(data);

   variable val_min = 1.0;
   variable val_max = 800;

   
   variable ywid = 8;
   variable pl = tikz_plot_new(13,ywid);
   pl.world(data[0].lo[0], data[0].hi[-1], 0.9*val_min, 1.1*val_max; loglog);

   variable cmap = qualifier("cmap","plasma");
   variable col_pal = get_color_palette(cmap, n);
   
   variable ii;
   _for ii(0, n-1){
      single_plot_data(pl, data[ii], sprintf("#%06x", col_pal[ii]), 0.7);
%      single_plot_data(pl, data[ii], "red", 0.7);
   }

   if (qualifier("noxtics",0)==1){
      pl.x1axis(;ticlabels=0);
  } else {
      pl.x1label(" Energy [keV]");
   }
   pl.y1label("Flux \; $[E \cdot F_E]$"R);
   
   variable p = tikz_plot_new(0.5, ywid*0.8);
   p.axis(; off);
   p.y2axis(; on);
   if (is_logdist(parvalues)) {
      p.world(0,1, min(parvalues), max(parvalues); ylog);
   } else {
      p.world(0,1, min(parvalues), max(parvalues));
   }
   p.y2label(strreplace(parname,"_", "\\_"));
   p.y2axis(; ticlabel_size="\footnotesize"R);
   p.plot_png(_reshape([0:255], [256, 1]); cmap=cmap, depth=200);
   p.plot([0,0,1,1,0], [0,1,1,0,0]; width=2, color="black", world0);
   
   return tikz_new_hbox_compound(pl, p, 0.5; center, interleave);
}

define plot_param(pname){

   variable data, parvalues;
   (data, parvalues) = eval_model_for_param(pname;; __qualifiers());

   return plot_data_struct(data, parvalues, pname);
}

variable ii, n = length(param_list);

variable pl_sum = Struct_Type[n];;

variable pl_list = {};
_for ii(0, n-1){
   pl_sum[ii] = plot_param(param_list[ii];log=is_log[ii], neval=neval_global);
   list_append(pl_list, pl_sum[ii] );
}

variable pl = tikz_new_vbox_compound(__push_list(pl_list), 1.2; interleave);
pl.render("fig_jedsad_relxill_comparison.pdf");
