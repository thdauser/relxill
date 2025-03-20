require("isisscripts");
require("tikz");


%%%
define get_eval_model(){
%%%
   variable emin = 0.3;
   variable emax = 70;
   variable nbins = 500;

   variable bin_lo, bin_hi;
   (bin_lo, bin_hi) = log_grid(emin, emax, nbins);

   variable emean = 0.5*(bin_lo + bin_hi);
   
   variable val = eval_fun_keV(bin_lo, bin_hi) / (bin_hi - bin_lo) * emean^2;

   return struct{elo=bin_lo, ehi=bin_hi, emean=emean, flux=val} ;
}


%%%
define get_parvalues(_vmin, _vmax, n){
%%%

   variable lo, hi;

   if (qualifier_exists("log") && qualifier("log",1)==1){
      (lo,hi) = log_grid(_vmin, _vmax, n-1);
   } else {
      (lo, hi) = linear_grid(_vmin, _vmax, n-1);
   }

   return [lo,hi[-1]];
}


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


%%%
define init_fitfun(_fitfun){
%%%
   fit_fun_default(_fitfun);
   set_par("*switch_reflfrac_boost",1);
}


%%%
define get_model_evaluations(_fitfun, parname, val_min, val_max){
%%%
   variable pl = tikz_plot_new();

   variable neval = qualifier("neval", 10);
   variable dat = Struct_Type[neval];
   
   variable parvalues = get_parvalues(val_min, val_max, neval ;;__qualifiers() );
   variable full_params = Array_Type[neval];

   init_fitfun(_fitfun);

   variable ii;
   _for ii(0, neval-1){
      set_par("*."+parname, parvalues[ii]);
      dat[ii] = get_eval_model();
      full_params[ii] = get_params();      
   }

   variable evalstruct =  struct_combine(struct{dat=dat}, struct{parvalues = parvalues} ) ;            
   variable auxinfo = struct{fitfun = _fitfun, parname=parname, full_params = full_params};
   
   return struct_combine(evalstruct, auxinfo);
}


%%%
define is_logdist(val){
%%%

   if (length(val) < 3 ){
      message(" *** error: need at least 3 model evaluation per parameter");
      return 0;
   }
   
   variable d1 = val[2] - val[1];
   variable d2 = val[1] - val[0];

   if ( abs(d1-d2) < 1e-4){
      return 0;
   } else {
      return 1;
   }
   
}


%%%
define plot_model_evaluations(dat){
%%%
   variable n = length(dat.dat);

   variable all_data = merge_struct_arrays(dat.dat);
   variable val_min = min(all_data.flux);
   variable val_max = max(all_data.flux);

   variable ywid = 8;
   variable pl = tikz_plot_new(13,ywid);
   pl.world(dat.dat[0].elo[0], dat.dat[0].ehi[-1], 0.9*val_min, 1.1*val_max; loglog);

   variable cmap = qualifier("cmap","plasma");
   variable col_pal = get_color_palette(cmap, n);
   
   variable ii;
   _for ii(0, n-1){
      pl.plot(dat.dat[ii].emean, dat.dat[ii].flux;
      opacity=0.4, color=sprintf("#%06x", col_pal[ii]));
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
   if (is_logdist(dat.parvalues)) {
      p.world(0,1, min(dat.parvalues), max(dat.parvalues); ylog);
   } else {
      p.world(0,1, min(dat.parvalues), max(dat.parvalues));
   }
   p.y2label(strreplace(dat.fitfun + "."+ dat.parname,"_", "\\_"));
   p.y2axis(; ticlabel_size="\footnotesize"R);
   p.plot_png(_reshape([0:255], [256, 1]); cmap=cmap, depth=200);
   p.plot([0,0,1,1,0], [0,1,1,0,0]; width=2, color="black", world0);

   
   return tikz_new_hbox_compound(pl, p, 0.5; center);
}

