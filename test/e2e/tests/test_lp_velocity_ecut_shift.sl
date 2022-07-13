require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";


define test_ener_shift(){

   variable lo,hi;
   variable val0, val1;
   (lo,hi) = log_grid(1,700,4000);

   set_par("relxilllpCp(1).switch_reflfrac_boost",1);
   
   set_par("*.Incl",20);
   set_par("*.kTe",40);
   set_par("*.gamma",2.0);
   set_par("*.switch_returnrad",0);

   variable veloc = [0.0,0.66];
   variable ii, n=length(veloc);
   
   % energies to calculate the flux in 
   variable elo = 100;
   variable ehi = 1000;
   variable ener_cutoff_p = Double_Type[n];
   variable ener_cutoff_r = Double_Type[n];
   if (n!=2){  % can only be 2 two as we compare beta=0 and beta>0
      return EXIT_FAILURE;
   }

   variable val_r = Array_Type[n];
   variable val_p = Array_Type[n];
   variable ener_fac = (0.5*(lo+hi))^2/(hi-lo);
   _for ii(0,n-1,1){
      set_par("*.beta",veloc[ii]);


      
      set_par("*.refl_frac",0,0,-1,10);
      val_p[ii] = eval_fun_keV(lo,hi)*ener_fac;
      set_par("*.refl_frac",-1,0,-1,10);
      val_r[ii] = eval_fun_keV(lo,hi)*ener_fac;

      val_p[ii] *= 1./val_p[ii][0];     
      val_r[ii] *= 1./val_r[ii][0];

      ener_cutoff_p[ii] = lo[wherefirst(val_p[ii] < 0.1)];
      ener_cutoff_r[ii] = lo[wherefirst(val_r[ii] < 0.1)];
      
      if (qualifier_exists("plot")){
	 title(" primary spectrum and reflection normalized to 1 at 1keV");
	 yrange(5e-3,50);
	 xlog;ylog;
	 line_style(1);color(ii+1);
	 ohplot(lo,hi,val_p[ii]);
	 color(ii+1);
	 oplot([ener_cutoff_p[ii],ener_cutoff_p[ii]],[1e-4,1e2]);

	 line_style(2);color(ii+1);
	 ohplot(lo,hi,val_r[ii]);
	 color(ii+1);
	 oplot([ener_cutoff_r[ii],ener_cutoff_r[ii]],[1e-4,1e2]);
	 line_style(1);

      }


      
      
   }

   
   if (not (ener_cutoff_p[0] < ener_cutoff_p[1]) ){
      msg_log += sprintf(" *** error: primary spectrum  kTe(beta=0)=%.3e < kTe(beta=%.2f)=%.3e not fulfilled",
		  ener_cutoff_p[0],veloc[1], ener_cutoff_p[1]);      
      return EXIT_FAILURE;
   }
   
   if (not (ener_cutoff_r[0] > ener_cutoff_r[1]) ){
      msg_log += sprintf(" *** error: primary spectrum  kTe(beta=0)=%.3e > kTe(beta=%.2f)=%.3e not fulfilled",
		  ener_cutoff_r[0],veloc[1], ener_cutoff_r[1]);
      return EXIT_FAILURE;
   }

   if (qualifier_exists("plot")){
      sleep(10);
      close_plot;
   }
   

   return EXIT_SUCCESS;
}

define runtest(ffs){ 
   
   msg_log += sprintf(" testing shift of Ecut for beta>0 \n");
   
   variable ff = "relxilllpCp";

   variable height = 30.0;
   variable ii,n = length(height);
   
   variable goodn;
   
   
   
   fit_fun_default(ff);
   set_par("*.h", height);
   
   return test_ener_shift();
   
}
