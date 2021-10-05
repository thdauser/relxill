require("load_test_setup.sl","Global");
% -*- mode: slang; mode: fold -*-

variable msg_log = "";

define runtest(ffs){ 
      
   variable status = EXIT_SUCCESS; 

   putenv("RELXILL_OUTPUT_FILES=1");
   
   fit_fun("relline_lp");
   set_par("*.h",3);
   
   () = eval_fun(1,2);
   variable rm,emis_mod_total;
   variable fname_output_relxill = "__relxillOutput_emisProfile.dat";
   (rm, emis_mod_total) = readcol(fname_output_relxill,1,2);
   
   set_par("*.switch_returnrad",0);   
   () = eval_fun(1,2);
   variable emis_mod_norrad;
   (rm, emis_mod_norrad) = readcol(fname_output_relxill,1,2);
 
   variable emis_mod_rrad = emis_mod_total - emis_mod_norrad;
   
   if (qualifier_exists("plot")){
      xlog;ylog;
      plot(rm, emis_mod_total);
      oplot(rm, emis_mod_norrad);
      oplot(rm, emis_mod_rrad);
      sleep(10);
   }
   
   
   if ( ( length (where(emis_mod_total <emis_mod_rrad)) !=0 )){
      message("   *** error: emissivity profile including return radiation does not have more flux in every bin than without");
      return EXIT_FAILURE;
   }

   variable rad_small = where(rm < 1.8);
   variable rad_large = where(rm > 100.0);

   %% direct radiation dominates for small radii
   if (  length(where(emis_mod_rrad[rad_small] > emis_mod_norrad[rad_small])) !=0 ){
      message("   *** error: emissivity profile not dominated at low energies by direct radiation ");
      return EXIT_FAILURE;
   }

   %% return radiation dominates for small radii
   if (  length(where(emis_mod_rrad[rad_large] < emis_mod_norrad[rad_large])) !=0 ){
      message("   *** error: emissivity profile not dominated at high energies by return radiation ");
      return EXIT_FAILURE;
   }

   
   putenv("RELXILL_OUTPUT_FILES=0");
   
   return status;
}

