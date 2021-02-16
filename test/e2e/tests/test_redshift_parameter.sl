require("load_test_setup.sl","Global");

variable msg_log = "";


define check_z_param(ff){ %{{{
   fit_fun(ff);
   
   variable zp = 0.1;
   
   variable lo, hi;
   (lo,hi) = log_grid(0.2,900,3000);
   
   set_par("*.z",0);
   variable val = eval_fun_keV(lo,hi);
      
   set_par("*.z",zp);
   variable val_z = eval_fun_keV(lo,hi);

   set_par("*.z",0.0);

   % now change grid
   variable zlo = lo * (1+zp) ;
   variable zhi = hi * (1+zp) ;
   
   variable val_z_corr = rebin(lo,hi,zlo,zhi,val_z);


   
   
   variable i_en;
   
   if (string_match(ff,"\.*line") == 0){
      i_en = where(1.0 < lo < 120.0);
   } else {
      i_en = where(1.0 < lo < 6.0);
   }
	
   lo = lo[i_en];
   hi = hi[i_en];
   zlo = zlo[i_en];
   zhi = zhi[i_en];
   val = val[i_en];
   val_z = val_z[i_en];
   val_z_corr = val_z_corr[i_en];
   
   variable ind = where(val_z_corr!=0);

   variable fac = (0.5*(lo+hi))/(hi-lo);
   
   if (not (qualifier_exists("nopl"))) {
      xlog;
      ylog;
      open_plot("/xw",2,1);
      hplot(lo,hi,val*fac);
      ohplot(lo,hi,val_z_corr*fac);
%      ohplot(zlo,zhi,val_z/(zhi-zlo)*(0.5*(zlo+zhi))^2);
      ohplot(lo,hi,val_z*fac);
      
      ylin;
      yrange(0.95,1.05);
      hplot(lo[ind],hi[ind],val[ind]/val_z[ind]);
      ohplot(lo[ind],hi[ind],val[ind]/val_z_corr[ind]);
   }
   
   variable gn = goodness(val,val_z_corr);
   msg_log += sprintf("    -> %s goodness value: %.3e\n",
	    ff, gn);
   return gn;
}
%}}}



define runtest(ffs){ 
   
   msg_log += sprintf(" testing REDSHIFT parameter\n");
   
   variable ff;   
   foreach ff(ffs){            
      if (check_z_param(ff;nopl) > goodness_lim){
	 msg_log += sprintf(" *** error: there seems to be a problem with the redshift in %s \n",ff);
	 return EXIT_FAILURE;
      }
   }

   return EXIT_SUCCESS;
}
%}}}


