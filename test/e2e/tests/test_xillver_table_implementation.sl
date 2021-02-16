require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";


define check_xilltab_implementation_single(ff,tabname){ %{{{
   
   
   %%%%  CURRENTLY NOT WORKED (KILLDE BY ISIS DUE TO MEMEORY (?) ISSUES)
   %% variable tablepath =  getenv("RELXILL_TABLE_PATH")+ "/";
   %% variable ff_tab = "tab";
   %% add_atable_model(tablepath+tabname,"tab");
   %% fit_fun(ff_tab);
   %% variable pars = get_params();
   %% set_params_xillver(pars);
   %% variable valr =  eval_fun_keV(lo0,hi0);
   %% valr *=  sum(val1) / sum(valr);

   variable XILLVER_REFDATA_DIR = "tests/refdata_xillverTable";
   
   %% doing this instead 
   variable refdat = fits_read_table(sprintf("%s/refdat_%s.fits",XILLVER_REFDATA_DIR,ff));
   variable lo0 = refdat.lo;
   variable hi0 = refdat.hi;
   variable valr = refdat.val;
   
   fit_fun_default(ff);
   set_par("*.refl_frac",-1.0);
   
   if (get_params("*.logxi")[0] != NULL){
      set_par("*.logxi",3.0); % we want to be on a grid point
   }
      
   if (string_matches(ff,"NS")!=NULL){
      set_par("*.kTbb",2.1);
   } else {
      set_par("*.gamma",2.1);
   }
      
   set_par("*.Incl",40);

   variable val1     =  eval_fun_keV(lo0,hi0);
   

   if ( goodness(val1,valr) > goodness_lim || qualifier_exists("pl")){
      xlog;ylog;
      hplot(lo0,hi0,val1);
      ohplot(lo0,hi0,valr);
      sleep(10);      
      hplot(lo0,hi0,val1/valr);
      sleep(10);      
   }
   
   return goodness(val1,valr);

}
%}}}


define runtest(ffs){ 
      
   
   msg_log += sprintf(" testing XILLVER models (compare to table models) \n");

   

   variable ff =     ["xillverCp", "xillver","xillverD","xillverNS","xillverCO"];
   variable ff_tab = ["xillver-comp.fits", "xillver-a-Ec5.fits", "xillverD-5.fits", "xillverNS.fits","xillverCO.fits"];

#ifdef STABLE
   ff =     ["xillverCp", "xillver","xillverD"];
   ff_tab = ["xillver-comp.fits", "xillver-a-Ec5.fits", "xillverD-5.fits"];   
#endif   
      
   variable ii,n = length(ff);
   
   _for ii(0,n-1,1){
      if (not (check_xilltab_implementation_single(ff[ii],ff_tab[ii] ) < goodness_lim)){
	 msg_log += sprintf(" *** error: MODEL %s does not agree with its table %s  \n",ff[ii],ff_tab[ii]);
	 return EXIT_FAILURE;
      }
      msg_log += sprintf("  comparing %s \t with table %s \t succesful \n", ff[ii], ff_tab[ii]);
   }
   
   return EXIT_SUCCESS;
}
%}}}

