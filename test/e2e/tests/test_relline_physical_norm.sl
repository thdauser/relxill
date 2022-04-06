require("load_test_setup.sl", "Global");

variable msg_log = "\n";

define runtest(ffs){
   
   variable ff;
   variable stdpar;
   foreach ff(ffs){
      
      %% does not make sense for xillver models
      if (string_match(ff,"rel")==0){
	 continue;
      }
      
      putenv("RELLINE_PHYSICAL_NORM");
      fit_fun(ff);
      variable val_notset = eval_fun_keV(0.05,2);
      
      stdpar = get_par("*.a");
      
      putenv("RELLINE_PHYSICAL_NORM=0");
      set_par("*.a",stdpar*0.99);
      () = eval_fun_keV(1,2);
      set_par("*.a",stdpar);
      variable val_set0 = eval_fun_keV(0.05,2);
       
      
      putenv("RELLINE_PHYSICAL_NORM=1");
      set_par("*.a",stdpar*0.98);
      () = eval_fun_keV(1,2);
      set_par("*.a",stdpar);
      variable val_set1 = eval_fun_keV(0.05,2);

      msg_log += sprintf(" *** %s (integ flux):  not set (%.5e) - ENV=0 (%.5e) - ENV=1 (%.5e)",
	       ff,sum(val_notset),sum(val_set0),sum(val_set1)) + "\n";
      
      variable prec = 1e-5;
      variable env0_changes_default = abs(sum(val_notset) - sum(val_set0)) > prec;
      variable env1_changes_default = abs(sum(val_notset) - sum(val_set1)) > prec;

      putenv("RELLINE_PHYSICAL_NORM");
            
      variable status = EXIT_SUCCESS;
      
      if (is_substr(ff, "relxilllp") ){	
	 %% should always have physical normalization, does not depend
	 %% on ENV
	 if (not ( env0_changes_default==0 && env1_changes_default==0 )) {
	    vmessage( "normalization of %s does not behave as expected", ff);
	    status = EXIT_FAILURE;
	 }
	 
      } else if( is_substr(ff, "relline") ){	 
	 %% by default normalized, but should change if ENV=1
	 if (not ( env0_changes_default==0 && env1_changes_default==1 ) ) {
	    status = EXIT_FAILURE;	    
	 }
	 
      } else if( is_substr(ff, "relxill") ){
	 %% by default normalized, but should change if ENV=1 [starting v1.4.1]
	 if ( not ( env0_changes_default==0 && env1_changes_default==1 )) {
	    status = EXIT_FAILURE;	    
	 }
	 
      } else {
	 vmessage(" ERROR: uncategorzied model %s", ff);
      }
      
      
      if (status==EXIT_FAILURE) {
	 vmessage(" *** error: normalization test failed!");
	 return EXIT_FAILURE;
      }
   }

   return EXIT_SUCCESS;
}


