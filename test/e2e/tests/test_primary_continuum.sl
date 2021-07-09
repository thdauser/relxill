require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";


define check_prim_cont_single(ff,ff_cont,assoc){ %{{{

   
   variable val0,val1;
   
   fit_fun_default(ff);
   if (isLpModel(ff)){
      set_par("*.boost",0.0,0,-10,10);      
   } else {
      set_par("*.refl_frac",0.0,0,-10,10);
   }
   val1 =  eval_fun_keV(lo0,hi0);

   variable z = 0.0; 
   if (string_match(ff,"lp")){
      z = grav_redshift_prim(get_par("*.a"),get_par("*.h"));
   }
   
   fit_fun_default(ff_cont);
   
   fit_fun(ff+"+"+ff_cont);
   variable key;
   foreach key(assoc_get_keys(assoc)){
      set_par("*."+assoc[key],get_par("*."+key));
   }
      
   fit_fun(ff_cont);
   
   % current definition for the high density models 
   if (string_match(ff,"D")){
      set_par("*.HighECut",300.0);
   }
   if (string_match(ff,"Cp")){
      set_par("nthComp*kT_bb",0.05);
      set_par("nthComp*inp_type",1.0);
      set_par("nthComp*Redshift",z);
   }

   val0 =  eval_fun_keV(lo0,hi0);
   val0 = val0/sum(val0)*sum(val1);
      
   if (not qualifier_exists("nopl")){
      xlog;ylog;
      hplot(lo0,hi0,val1);
      ohplot(lo0,hi0,val0);
   }
   
   variable gn = goodness(val1,val0/sum(val0)*sum(val1));
   msg_log += sprintf("   -> primary continum %s for reflection model %s  [gn=%.2e]\n",ff_cont,ff,gn);

   return gn; 
}
%}}}



define runtest(ffs){ 

   msg_log += sprintf("testing PRIMARY CONTINUUM\n",counter);

   
   variable assoc = Assoc_Type[String_Type];
   assoc["gamma"] = "PhoIndex";
   assoc["Ecut"]  = "HighECut";
   variable assocD = Assoc_Type[String_Type];   
   assocD["gamma"] = "PhoIndex";
   variable assocCp = Assoc_Type[String_Type];
   assocCp["gamma"] = "Gamma";
   assocCp["kTe"] = "kT_e";
   variable assocBB = Assoc_Type[String_Type];
   assocBB["kTbb"] = "kT";

   variable ff =      ["relxill","relxilllp","relxillD","relxilllpD","relxillCp","relxilllpCp","xillver","xillverD","xillverCp","relxillNS","xillverNS"];
   variable ff_cont = ["cutoffpl","cutoffpl","cutoffpl","cutoffpl","nthComp","nthComp","cutoffpl","cutoffpl","nthComp","bbody","bbody"];
   variable arr_assoc=[assoc,assoc,assocD,assocD,assocCp,assocCp,assoc,assocD,assocCp,assocBB,assocBB];
   
   
   variable ii,n = length(ff);
   
   _for ii(0,n-1,1){

      if ( model_exists(ff[ii], ffs)) {      
	 if (not (check_prim_cont_single(ff[ii],ff_cont[ii],arr_assoc[ii]; nopl) < goodness_lim)){
	    msg_log += sprintf(" *** error: there seems to be a problem with the PRIMARY CONTINUUM in  MODEL %s \n",ff[ii]);
	 return EXIT_FAILURE;
	 }
      } 
	     
   }
   
   return EXIT_SUCCESS;
}
