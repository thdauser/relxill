require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";


define is_rrad_implemented(ff){
   
   fit_fun_default(ff);
   
   variable status = 1;
   
   if (get_params("*.switch_returnrad")[0]==NULL){
      vmessage("    model %s does not have returning radiation implemented", ff);
      status = 0;
      
      % require that the model with returning radiation has a larger
      % flux (for relxill models)
   } else if (is_relxill_model(ff)){
	 
      variable val_rrad = eval_fun_keV(lo0, hi0);
      set_par("*.switch_returnrad", 0);
      variable val_no_rrad = eval_fun_keV(lo0, hi0);
      
      if ( not  (sum(val_rrad) / sum(val_no_rrad) > (1 + 1e-6) ) ) {
	 
	 vmessage("  *** error:  model %s shows not a larger flux of returning radiation is included", ff);
	 vmessage("        (rrad:%e  no-rrad:%e)      ", sum(val_rrad) , sum(val_no_rrad) );
	 
	 status = 0;
      }
      
      
   }
	      
   return status;
	  
}

define runtest(ffs){ 
      
   variable status = EXIT_SUCCESS; 
   
   %  test that for the given models RETURN RADIATION is implemented
   variable fit_fun_rrad = 
     ["relline_lp", "relxilllp", "relxilllpCp"];
   
   
   variable ii, n = length(fit_fun_rrad);
   variable implem = Int_Type[n];
     
   _for ii(0, n-1){
      implem[ii] = is_rrad_implemented(fit_fun_rrad[ii]);
      
   }
   
   if (sum(implem) != n){
      status = EXIT_FAILURE;
   }
   
   return status;
}

