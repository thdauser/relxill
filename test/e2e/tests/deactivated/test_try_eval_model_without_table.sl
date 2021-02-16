
define runtest(ffs){ %{{{


   variable tab_env = getenv("RELXILL_TABLE_PATH");
   putenv("RELXILL_TABLE_PATH");

   msg_log = "\n";
   
   variable ff;
   variable val;
   foreach ff(ffs){
      msg_log += sprintf(" trying to load %s without table ",ff);
      fit_fun(ff);
      try{
	 val = eval_fun_keV(1,2);
      }catch AnyError;
      if (not ( val >= 0 )){
	 vmessage(" *** error: simple test withOUT tables for %s failed!",ff);
	 return EXIT_FAILURE;
      }
   }

   putenv(sprintf("RELXILL_TABLE_PATH=%s",tab_env));

   return EXIT_SUCCESS;
}
%}}}
