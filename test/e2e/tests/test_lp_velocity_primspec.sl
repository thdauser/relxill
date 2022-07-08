require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";

variable lo,hi;
variable val0, val1;
(lo,hi) = log_grid(1,700,4000);
variable en = 0.5*(lo+hi);

define gam(bet) {
      return 1./sqrt(1.-bet^2);
}

define doppler(theta,bet) {
      return 1./(gam(bet)*( 1. + bet*cos(theta) )  );
}



define test_veloc_boost(ff,veloc){

  
   fit_fun_default(ff);
   set_par("relxilllpCp(1).switch_reflfrac_boost",1);   
   set_par("*.Incl",20);
   set_par("*.kTe",40);
   set_par("*.gamma",2.0);
   
   set_par("*.h", 80);

   set_par("relxilllpCp(1).switch_reflfrac_boost",1);

   
   set_par("*.Incl",30);
   set_par("*.kTe",40);
   set_par("*.gamma",2.0);
   
   variable ii, n=length(veloc);
   variable val_r = Array_Type[n];
   variable val_p = Array_Type[n];
   variable ratio_sum_predicted_boost = Double_Type[n];
   _for ii(0,n-1,1){
      set_par("*.beta",veloc[ii]);

      
      set_par("*.refl_frac",0);
      val_p[ii] = eval_fun_keV(lo,hi)/(hi-lo);

      
      variable ener_shift_0 = 1/(1+kerr_lp_redshift(get_par("*.h")[0], get_par("*.a")[0]));
      variable doppler_factor =  doppler(get_par("*.Incl")[0]*PI/180., -veloc[ii]);
      
      variable boost = (ener_shift_0*doppler_factor) ^ get_par("*.gamma")[0] * doppler_factor ^ 2;
      variable sum_norm = sum(val_p[ii][where( 1<= lo <= 10)]);

      ratio_sum_predicted_boost[ii] = sum_norm/boost;
	
      msg_log += sprintf("beta:%.2f --  sum in 1-10keV: %.4e (ener shift: %.4e, doppler_factor: %.4e, boost PL:%.4e) norm/boost:%.3e",
	       veloc[ii], sum_norm, ener_shift_0, doppler_factor, boost, sum_norm/boost
	      );
      
      
%      yrange(5e-6,500);
%      xlog;ylog;
%      line_style(1);color(ii+1);
%      ohplot(lo,hi,val_p[ii]);
   }

   
   %% the calculated boost should be rather close to the change in flux 
   variable relative_diff = abs(ratio_sum_predicted_boost/mean(ratio_sum_predicted_boost)-1);
   if ( length(where(relative_diff > 0.05)) == 0) {
      return EXIT_SUCCESS;
   } else {
      return EXIT_FAILURE;
   }
   
   
}


define plot_height(height){

%%   putenv("RELXILL_PRINT_DETAILS=1");


   
   variable ii, n=length(height);
   variable val_r = Array_Type[n];
   variable val_p = Array_Type[n];
   _for ii(0,n-1,1){
      set_par("*.h",height[ii]);

      
      set_par("*.refl_frac",0);
%      val_p[ii] = eval_fun_keV(lo,hi)/(hi-lo)*(en)^2;
      val_p[ii] = eval_fun_keV(lo,hi);

      variable ener_shift = 1/(1+kerr_lp_redshift(height[ii], get_par("*.a")[0]));
      variable boost = ener_shift ^ get_par("*.gamma")[0];
      
      msg_log += sprintf("height:%.2f --  sum in 1-10keV: %.4e (ener shift: %.4e, boost PL:%.4e)",
			 height[ii], sum(val_p[ii][where( 1<= lo <= 10)]), ener_shift, boost);
      
      variable norm_fac = 1.0; % 1./val_p[0][0];;
      
      val_p[ii] *= norm_fac;
      
      yrange(5e-6,500);
      xlog;ylog;
      line_style(1);color(ii+1);
      ohplot(lo,hi,val_p[ii]/(hi-lo)*(en)^2);
   }

}


define runtest(ffs){ 
   
   msg_log += sprintf(" testing LAMP POST VELOCITY\n");
   
   variable ff = "relxilllpCp";

   variable height = [30.0,5.0,3.0];
   variable veloc = [0.0,0.01,0.1,0.2,0.3];
   variable ii,n = length(height);
   
   
   
   if (test_veloc_boost(ff,veloc) != EXIT_SUCCESS)
     return EXIT_FAILURE;
   

   
   _for ii(0,n-1,1){

%      vmessage(" height = %.2f ", height[ii] );
%      set_par("*.h", height[ii]);
%      plot_veloc(veloc);

%      set_par("*.beta",0.0);
%      plot_height(height);

      
   }

   
   return EXIT_SUCCESS;
}
