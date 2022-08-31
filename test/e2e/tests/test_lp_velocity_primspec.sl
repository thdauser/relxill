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


define set_model(ff){

   fit_fun_default(ff);
   set_par("relxilllpCp(1).switch_reflfrac_boost",1);   
   set_par("*.Incl",20);
   set_par("*.kTe",40);
   set_par("*.gamma",1.7);
   
   set_par("*.h", 50);

   set_par("relxilllpCp(1).switch_reflfrac_boost",1);

   

   
}


define calc_boost_beta_primary(beta){

   variable doppler_factor =  doppler(get_par("*.Incl")[0]*PI/180., -beta);

   variable gam = get_par("*.gamma")[0];
   variable boost = (doppler_factor) ^ gam * doppler_factor ^ 2;
   
   msg_log += sprintf("beta:%.2f --  (doppler_factor: %.2f, boost PL:%.2f) ",
		      beta, doppler_factor, boost    );


   return boost;
}



define calc_boost_beta_reflected(beta){

   variable gam = get_par("*.gamma")[0];
   variable boost = ( 1-(beta+1)^(-1-gam) ) / (beta*(1+gam)) * (1-beta^2)^(1+gam/2);
   
   msg_log += sprintf("beta:%.2f -- reflected boost PL:%.4e ", beta, boost);

   return boost;
}


define flux_in_band(beta, refl_frac){

   variable emin = 1.0;
   variable emax = 10.0;
   
   set_par("*.beta",beta);
   set_par("*.refl_frac",refl_frac,0,-1,10);
   
   variable val = eval_fun_keV(lo,hi)/(hi-lo);
   variable sum_norm = sum(val[where( emin<= lo <= emax)]);

   
   if (qualifier_exists("plot")){
      yrange(5e-2,3500);
      
      xlog;ylog;
      line_style(nint(refl_frac)+2);
      ohplot(lo,hi,val*(0.5*(lo+hi)^2));
   }

   return sum_norm;
}


%%%%%%%%%%%
define test_veloc_boost(ff,veloc){
%%%%%%%%%%%

   set_model(ff);

   variable refl_frac = 0.0;
   if (qualifier_exists("reflected"))
     refl_frac = -1.0;
      
   variable sum_norm_0 = flux_in_band(0.0, refl_frac);

   color(1);
   
   variable ii, n=length(veloc);
   variable ratio_sum_predicted_boost = Double_Type[n];
   _for ii(0,n-1,1){

      variable boost =
	(qualifier_exists("reflected"))
	? calc_boost_beta_reflected(veloc[ii])
	: calc_boost_beta_primary(veloc[ii]);
      variable sum_norm = flux_in_band(veloc[ii], refl_frac) / sum_norm_0;
      
      ratio_sum_predicted_boost[ii] = sum_norm/boost;
      
      msg_log += sprintf("          sum in 1-10keV: %.4e -> norm/boost:%.2f\n", sum_norm, ratio_sum_predicted_boost[ii]);            
   }

   
   %% the calculated boost should be rather close to the change in flux 
   variable reference_value = 1.0;
   variable relative_diff = abs(ratio_sum_predicted_boost - reference_value);
   if ( length(where(relative_diff > 0.1)) == 0) {
      return EXIT_SUCCESS;
   } else {
      vmessage(" *** error *** ");
      return EXIT_FAILURE;
   }
   
   
}



%%%%%%%%%%%%%%%%%
define runtest(ffs){ 
%%%%%%%%%%%%%%%%%
   
   
   variable ff = "relxilllpCp";

   variable veloc = [0.01,0.1,0.2,0.3,0.66];
   
   
   msg_log += sprintf(" testing LAMP POST VELOCITY BOOST for PRIMARY SOURCE\n");
   
   if (test_veloc_boost(ff,veloc;primary) != EXIT_SUCCESS)
     return EXIT_FAILURE;
   
   msg_log += sprintf("\n testing LAMP POST VELOCITY BOOST for REFLECTED SPECTRUM (larger uncertainty)\n");

   if (test_veloc_boost(ff,veloc;reflected) != EXIT_SUCCESS)
     return EXIT_FAILURE;
   
   return EXIT_SUCCESS;
}
