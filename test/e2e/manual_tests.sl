define do_mc_testing(){ %{{{

   
   counter++;
   vmessage("\n### %i ### random MC Parameter tests ### ", counter);

   variable n_mc = 100;
   variable stat;

   variable lo,hi;
   (lo,hi) = log_grid(0.1,500,500);
   variable fakeval = lo^-(1.9);
   variable fake_dat = struct{
         bin_lo=lo,
         bin_hi=hi,
         value=fakeval,
         err= 0.1*fakeval
   };
   
   variable iDat = define_counts(_A(fake_dat));

   variable save_fit_verbose = Fit_Verbose;
   Fit_Verbose=0;
   
   
   variable ffs = ALL_FF;
   variable ff;
   variable val, dt;
   foreach ff(ffs){

      fit_fun(ff);
      freeze("*.norm");

      tic;
      stat = fit_search(n_mc, &eval_counts;  serial, dir="/tmp/fit_search");
      dt = toc / n_mc * 1e3;
      if ( stat == NULL ){
	 vmessage(" *** error: MC Parameter test for %s failed!",ff);
	 return EXIT_FAILURE;
      }
      vmessage("   -> tested %s with %i evaluations   \t [ <time> ~ %.0fmsec ]",ff,n_mc,dt);

   }
   
   Fit_Verbose = save_fit_verbose;
   delete_data(iDat);
   
   return EXIT_SUCCESS;
}
%}}}

define print_refl_frac(){ %{{{
   
   counter++;
   vmessage("\n### %i ### print REFLECTION FRACTION information: ###",counter);
   
   variable ff = ["relxilllp","relxilllpD","relxilllpCp"];
   
   variable ii,n = length(ff);
   
   variable goodn;
   _for ii(0,n-1,1){
      fit_fun(ff[ii]);
      vmessage(" %s : ",ff[ii]);
      set_par("*.fixReflFrac",2);
      () = eval_fun_keV(1,2);
      set_par("*.fixReflFrac",0);
      message("\n");
   }
   
   return EXIT_SUCCESS;
}
%}}}

print_refl_frac();
do_mc_testing();
