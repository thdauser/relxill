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


define getRgridFromHeader(fname){   
   return fits_read_key(fname,"rlo","rhi");   
}

define testSingleXillverPrimaryNormalization(inp_struct, izones){
    
    variable pl = xfig_plot_new(18,12);
    
    pl.world(0.1,50,1e-9,1;loglog);
    
    pl.title("\noindent Need new \\"R+
    +"description \\"R);

%    variable ii, n = length(dat);  %% need to change that
%    variable xillNorm = Double_Type[n];
    
   print(inp_struct);
   
    % variable col_ind = [0,1,1,3,3,6];
    % variable line_arr = [0,0,1,0,1,0];
    % variable line_width=[3,1,1,1,1,3];
    
    % _for ii(0,n-1){
       
    %   variable en = 0.5*(dat[ii].bin_lo+dat[ii].bin_hi);
    % 	pl.plot(en, dat[ii].flux;width=line_width[ii], color=CB_COLOR_SCHEME_NB[col_ind[ii]], line=line_arr[ii]);
	
    % 	xillNorm[ii] = get_xillver_norm(en, dat[ii].flux);
	
    % 	vmessage(" norm flux wrt xillver %s: %e", fname[ii], xillNorm[ii]);
    % }
    
    % pl.xlabel("Energy [keV]");
    % pl.ylabel(" Counts / Bins ");
    
    % variable rlo, rhi;
    % pl.add_object(
    % xfig_new_text(sprintf("$T_\mathrm{in}=%.1f$\,keV;  Radial Zone #%i"R,Tin,izone)),
    % 0.96,0.96,0.5,0.5;world0
    % );
        
    variable retVal = EXIT_SUCCESS;
    
   %  if ( abs(xillNorm[0]-xillNorm[2]) > 1e-6){
   % 	retVal = EXIT_SUCCESS;
   % 	vmessage(" *** warning : incident and primary xillver spectrum do not have the same normaliztion [%e,%e]",
   % 	xillNorm[0],xillNorm[2]);
   % }
   
   % if ( xillNorm[1]>=xillNorm[2]){
   %     retVal = EXIT_SUCCESS;
   %     vmessage(" *** warning : reflected xillver norm is larger than primary xillver norm [%e,%e]",
   %     xillNorm[1],xillNorm[2]);
   % }      
   
   return (retVal,pl);
}



define createOutputFiles(){ 
   load_xspec_local_models("build/");
   fit_fun(global_fit_fun);
   
   putenv("RELXILL_WRITE_FILES=1");

   () = system("rm -f *");
   
   message("  Creating Output Files: ");

   list_par;
   
   set_par("*.a",spin);
   set_par("*.kTbb", Tin);
   set_par("*.Rout",400);
   set_par("*.boost",-1);
   set_par("*.Incl",40);
   () = eval_fun(1,2);

   set_par("*.boost",-1);
   () = eval_fun(1,2);
   
}
%}}}



define loadOutputFilesRframe() {
   variable rframe = getFilenameStruct().rframe;
   return struct  {
      xill_refl = get_2d_data(rframe.xillRefl,spin;diskarea=1),
      xill_prim = get_2d_data(rframe.xillPrim,spin;diskarea=1),
      spec_ret  = get_2d_data(rframe.specRet,spin;diskarea=1),
      spec_prim = get_2d_data(rframe.specPri,spin;diskarea=1)
   };
}


define testXillverPrimaryNormalization(){
   
   variable izones = [10,20,30,40];

   createOutputFiles();
   variable dat = loadOutputFilesRframe();

   variable retVal, pl;
   (retVal, pl) = testSingleXillverPrimaryNormalization(dat,izones);

   pl.render(dir_plots+"testXillverPrimaryNormalization.pdf";Verbose=-1);
   
   
   if (sum(retVal) == EXIT_SUCCESS){
      return EXIT_SUCCESS;
   } else {
      return EXIT_FAILURE;
   }

}



if (length(__argv)>0){
   () = testXillverPrimaryNormalization();
}
