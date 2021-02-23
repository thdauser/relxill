require("isisscripts");

load_xspec_local_models("build/");

variable lo,hi;
(lo,hi) = log_grid(0.5,9.5,500);

variable n = 15;
variable a = [0.8:0.99:#n];

variable pl = xfig_plot_new();
pl.world(0.5,9.5,-.3,.3);

fit_fun("rellinelpRet");

variable ii;
variable val = Array_Type[n];
_for ii(0,n-1){
   val[ii] = eval_fun_keV(lo,hi)/(hi-lo)*lo;
 
   variable rat = val[ii]-val[0];
   variable ind = where(val[0]>0);
   
   set_par("*.a",a[ii]);
   pl.plot(lo[ind],rat[ind];color=ii+1);
}

pl.render("plots/rellinelpRet_spin.pdf");