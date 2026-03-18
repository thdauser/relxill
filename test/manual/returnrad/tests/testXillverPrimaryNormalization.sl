require("isisscripts");

require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

%%%%%%%%%%%%%%%%%%
variable spin = 0.998;
variable Tin = 2.0;
variable logN = 18.0;
variable logxi = 2.0;

variable global_fit_fun = "relxillBBxill";
%%%%%%%%%%%%%%%%%%

require("create_output_files.sl");

variable ff_plots = "relxillBBxill";
variable prefix = ff_plots+"_";

%% =1= set parameters
variable spin = 0.998;  %% currently only done for this value
variable Tin = 2.0;
variable logxi = 2.0;
variable logn = 18.0;



define testSingleXillverPrimaryNormalization(inp_struct, izone){
    
    variable pl = xfig_plot_new(18,12);
    
    pl.world(0.1,50,1e-9,1;loglog);
    
    pl.title("\noindent Need new \\"R+
    +"description \\"R);

%    variable ii, n = length(dat);  %% need to change that
%    variable xillNorm = Double_Type[n];

   variable en = inp_struct.xill_refl.dat.ener[0,*];
   
   variable spec_prim = inp_struct.spec_ret.dat.spec[izone,*];
   variable spec_prim_xill = inp_struct.xill_prim.dat.spec[izone,*];
   variable spec_refl = inp_struct.xill_refl.dat.spec[izone,*];
   
   variable xillNorm_prim = get_xillver_norm(en, spec_prim);
   variable xillNorm_prim_xill = get_xillver_norm(en, spec_prim);
   variable xillNorm_refl = get_xillver_norm(en, spec_refl);

   if (qualifier_exists("verbose")){
      vmessage(" prim: %.3e, prim_xill: %.3e, refl: %.3e, ratio:%.3e", 
      xillNorm_prim, xillNorm_prim_xill, xillNorm_refl,xillNorm_prim_xill/xillNorm_refl);
   }


   variable retVal = EXIT_SUCCESS;
   if ( abs(xillNorm_prim-xillNorm_prim_xill) > 1e-6){
      retVal *= EXIT_FAILURE;  % as EXIT_FAILURE==0, this will make it fail overall
      vmessage(" *** test failure : incident and primary xillver spectrum do not have the same normaliztion [%e,%e]",xillNorm_prim-xillNorm_prim_xill);
   }
   
   if ( xillNorm_refl>=xillNorm_prim){
      retVal *= EXIT_FAILURE;
      vmessage(" *** test failure : reflected xillver norm is larger than primary xillver norm [%e,%e]",
      xillNorm_refl, xillNorm_prim);
   }      
   
   return retVal;
}


define loadOutputFilesRframe() {
   variable rframe = getFilenameStruct().rframe;
   return struct  {
      xill_prim = get_2d_data(rframe.xillPrim,spin;diskarea=1),
      xill_refl = get_2d_data(rframe.xillRefl,spin;diskarea=1),
      spec_ret  = get_2d_data(rframe.specRet,spin;diskarea=1),
      spec_prim = get_2d_data(rframe.specPri,spin;diskarea=1)
   };
}





define testXillverPrimaryNormalization(){

   
   
   %% =2= create the output files, pre-defined set of parameters
   variable params = get_param_struct(spin, Tin, logxi, logn);
   createOutputFiles(ff_plots,params);

   
   variable dat = loadOutputFilesRframe();

   variable izones = [0,1,5,10,30];
   variable ii, n = length(izones);
   variable retVal = Int_Type[n];

   _for ii(0, n-1){
      retVal[ii] = testSingleXillverPrimaryNormalization(dat,izones[ii];verbose);
   }
   
   if (sum(retVal) == EXIT_FAILURE){
      return EXIT_FAILURE;
   } else {
      return EXIT_SUCCESS;
   }
}



if (length(__argv)>1){
   () = testXillverPrimaryNormalization();
}
