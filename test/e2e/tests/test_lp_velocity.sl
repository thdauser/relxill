require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";


define plot_veloc(veloc){

   variable lo,hi;
   variable val0, val1;
   (lo,hi) = log_grid(1,700,4000);

   set_par("relxilllpCp(1).switch_reflfrac_boost",1);
   
   set_par("*.Incl",20);
   set_par("*.kTe",40);
   set_par("*.gamma",1.7);
   set_par("*.switch_returnrad",0);

   
   variable ii, n=length(veloc);
   variable val_r = Array_Type[n];
   variable val_p = Array_Type[n];
   _for ii(0,n-1,1){
      set_par("*.beta",veloc[ii]);

      set_par("*.refl_frac",0);
      val_p[ii] = eval_fun_keV(lo,hi)/(hi-lo)*(0.5*(lo+hi))^2;
      set_par("*.refl_frac",-1,0,-1,10);
      val_r[ii] = eval_fun_keV(lo,hi)/(hi-lo)*(0.5*(lo+hi))^2;

      variable norm_fac = 1./val_p[ii][0];
      variable norm_fac_r = 1./val_p[ii][0];


      if (qualifier_exists("plot")){
	 yrange(5e-4,5);
	 xlog;ylog;
	 line_style(1);color(ii+1);
	 ohplot(lo,hi,val_p[ii]*norm_fac);
	 line_style(2);color(ii+1);
	 ohplot(lo,hi,val_r[ii]*norm_fac_r);
	 line_style(1);
      }
   }

}

define runtest(ffs){ 
   
   msg_log += sprintf(" testing LAMP POST VELOCITY\n");
   
   variable ff = "relxilllpCp";

   variable height = [3.0];
   variable veloc = [0.0,0.44,0.88];
   variable ii,n = length(height);
   
   variable goodn;
   
   
   
   fit_fun_default(ff);
   _for ii(0,n-1,1){

%       vmessage(" height = %.2f ", height[ii] );
      set_par("*.h", height[ii]);
      plot_veloc(veloc);      
   }

   close_plot;
   
   return EXIT_SUCCESS;
}
