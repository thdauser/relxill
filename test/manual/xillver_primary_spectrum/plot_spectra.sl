require("isisscripts");
() = evalfile("load_model.sl");

variable logxi = [0:4.7:#6];
variable incl = [10:80:#6];

variable ff = "xillver";
variable cols = CB_COLOR_SCHEME_NB;

define set_xill_param(logxi, incl, refl_frac){
   fit_fun(ff);
   set_par("*.logxi",logxi);
   set_par("*.Incl",incl);
   set_par("*.refl_frac",refl_frac);
   
   
}


variable emin = 0.1;
variable emax = 1000;

define eval_xill_model(logxi, incl){
 
   variable lo, hi;
   (lo, hi) = log_grid(emin,emax,3000);
     
   
   set_xill_param(logxi, incl, -1);
   variable val_r = eval_fun_keV(lo,hi)/(hi-lo)*lo*lo;
   
   set_xill_param(logxi, incl, 0);
   variable val_p = eval_fun_keV(lo,hi)/(hi-lo)*lo*lo;
      
   return struct{lo=lo, hi=hi, prim=val_p, refl=val_r};
}

define get_std_plot(){
   variable pl = xfig_plot_new();
   pl.world(emin,emax, 1e-2,8e2; loglog);
   
   pl.xlabel("Energy [keV]");
   pl.ylabel("$\nu F_\nu$ Flux"R);
   
   return pl;
}


define plot_xill_logxi(logxi_arr, incl){
   
 
   variable ii, n = length(logxi_arr);
   
   variable pl = get_std_plot();
   
   _for ii(0,n-1){
      variable dat = eval_xill_model(logxi_arr[ii], incl);
      
      pl.plot(0.5*(dat.lo+dat.hi), dat.refl; width=1, color=cols[ii]);
      pl.plot(0.5*(dat.lo+dat.hi), dat.prim; width=2, color="black");
      
   }
   
   pl.add_object(xfig_new_text(sprintf("$i = %.0f^\circ$"R,incl)),
		 0.5,0.98,0,0.5;world0);
   
   return pl;
}

define get_label(vals, str){
   variable ii, n = length(vals);
   variable lab =String_Type[n];
   
   _for ii(0,n-1){
      lab[ii] = sprintf(str,vals[ii]);
   }
   return lab;
}

define plot_xill_spectra(logxi_arr, incl_arr){
   
   variable ii, n = length(incl_arr);
   variable pl = Struct_Type[n];
   
   _for ii(0,n-1){
      pl[ii] = plot_xill_logxi(logxi_arr,incl_arr[ii]);
   }
   

   variable lab = get_label(logxi_arr,"$\xi = %.1f$"R);
   pl[0].add_object(xfig_new_legend(lab,cols,0,1,0.5;labelsize="footnotesize"),
		   0.98,0.02,0.5,-0.5;world0);
   
   return xfig_multiplot(pl; cols=2);
}

plot_xill_spectra(logxi, incl).render("plot_all_xillspec.pdf");