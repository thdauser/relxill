require("load_test_setup.sl","Global");

variable msg_log = "";

define runtest(ffs){ 
      
   variable status = EXIT_SUCCESS; 

   variable ff;
   variable val,val2,val3;
   variable stdpar;
   foreach ff(ffs){
      
      if (string_match(ff,"rel")==0){
	 continue;
      }
      
      fit_fun(ff);
      variable val_notset = eval_fun_keV(1,2);
      
      variable prec = 1e-5;
      variable different_from_unity = abs(sum(val_notset) - 1.0 ) > prec;

      msg_log += sprintf(" - %s, different from unity: %.3e\n", ff, sum(val_notset));


      
      if (is_substr(ff, "relxill") ){	
	 %% should always have physical normalization, does not depend
	 %% on ENV
	 if ( different_from_unity==0 ) {
	    status = EXIT_FAILURE;
	 }
	 
      } else if( is_substr(ff, "relline") ){	 
	 %% by default normalized, but should change if ENV=1
	 if ( different_from_unity==1  ) {
	    status = EXIT_FAILURE;	    
	 }
	 
	 
      } else {
	 vmessage(" ERROR: uncategorzied model %s", ff);
      }

            
      if ( status==EXIT_FAILURE ){
	 vmessage(" *** error: normalization test failed!");
	 return status;
      }
   }
   
   return status;
}
%}}}
