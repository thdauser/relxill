require("load_test_setup.sl","Global");

variable msg_log = "";

define runtest(ffs){ %{{{

      
   variable ff;
   foreach ff(ffs){
      fit_fun(ff);

      msg_log += ff + "\n" ;

      variable val = eval_fun_keV(1,2);
      if (not ( val > 0 )){
	 vmessage(" *** error: simple test for %s failed! (val=%e)",ff,val);
	 return EXIT_FAILURE;
      }
   }
   
   return EXIT_SUCCESS;
}
%}}}
