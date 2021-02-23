require("isisscripts");

load_xspec_local_models("build/");

variable lo,hi;
(lo,hi) = log_grid(0.5,500,3000);

variable a = [0.85,0.9,0.99,0.998];
variable n = length(a);

variable pl = xfig_plot_new();
pl.world(0.5,300,1,300;loglog);

variable pl2 = xfig_plot_new();
pl2.world(0.5,300,0.37,1.2;loglog);

variable plr = xfig_plot_new();
%plr.world(0.5,500,0.02,2;loglog);
plr.world(0.5,300,0.79,1.21;xlog);

plr.xlabel("Energy [keV]");
plr.ylabel("Ratio [normalized]");

fit_fun("relxilllpRet");

variable ii;
variable val = Array_Type[n];
variable val0 = Array_Type[n];
_for ii(0,n-1){
   set_par("*.a",a[ii]);
   set_par("*.h",6,0,2,100);

   set_par("*.return_rad",0);
   set_par("*.refl_frac",-1,0,-1,10);
   set_par("*.fixReflFrac",3);
   val0[ii] = eval_fun_keV(lo,hi)/(hi-lo)*lo*lo;
   pl.plot(lo, val0[ii];color=ii+1,line=3);

   set_par("*.return_rad",1);
   val[ii] = eval_fun_keV(lo,hi)/(hi-lo)*lo*lo;    
   pl.plot(lo, val[ii];color=ii+1);


   
   pl2.plot(lo,val[ii]/val[ii][0]; color=ii+1);
   pl2.plot(lo,val0[ii]/val0[ii][0]; color=ii+1, line=1);

   variable rat = val[ii]/val0[ii];
   plr.plot(lo,rat/rat[0]  ; color=ii+1);
%   plr.plot(lo,val0[ii]/max(val0[ii])*2; color=ii+1, line=1);

   
   set_par("*.refl_frac",0);
   set_par("*.fixReflFrac",0);
   val[ii] = eval_fun_keV(lo,hi)/(hi-lo)*lo*lo;    
   pl.plot(lo, val[ii];color=ii+1,line=1);

   
   
}

xfig_multiplot(pl, pl2, plr).render("plots/relxilllpRet_spin.pdf");
