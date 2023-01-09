require("xspec");
load_xspec_local_models("build/");
require("isisscripts");

variable lo,hi;
(lo,hi) = log_grid(0.1,1000,4000);


define plotmodels(ff,pname){


   fit_fun(ff);
   list_par;

   variable efac = 1./(hi-lo)*(0.5*(lo+hi))^2;

   variable ecut = [100,40,10,400];
   variable ii, n=length(ecut);
   variable val_cont0 = Array_Type[n];

   xlog;ylog;
   yrange(1e-3,3e2);

   _for ii(0,n-1){      
      set_par(pname,ecut[ii]);
      
%      set_par("*.refl_frac",0);
      color(ii+1); line_style(1);

      val_cont0[ii] = eval_fun_keV(lo,hi)*efac;
      ohplot(lo,hi, val_cont0[ii]);

      color(ii+1); line_style(2);
      variable val_cont = eval_fun_keV(lo,hi)*efac;
      variable renorm = val_cont0[0][0]/val_cont[0];
      ohplot(lo,hi, val_cont*renorm);
     
      

%      set_par("*.refl_frac",-1);
%      color(ii+1); line_style(2);
%      ohplot(lo,hi,eval_fun_keV(lo,hi)*efac);
   }
     
   sleep(100);

}


   
%variable ff =      "xillverCp";
variable ff =      "nthcomp";

plotmodels(ff,"*.kT_e");

   
   
