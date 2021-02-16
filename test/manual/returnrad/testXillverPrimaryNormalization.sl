require("isisscripts");

require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");

define getRgridFromHeader(fname){
   
   return fits_read_key(fname,"rlo","rhi");
   
}

define testSingleXillverPrimaryNormalization(izone){
   
   variable noTshiftString  = "Tshift1.00";
   variable hasTshiftString = "Tshift1.40";
   
   variable dir = "build/";
   variable fname = dir+ "testrr-spec-rframe-bbody-" + sprintf("izone%02i",izone) + 
     [
      "-$noTshiftString-incident.fits"$,
      "-$noTshiftString-reflect.fits"$,
      "-$noTshiftString-reflectPrim.fits"$,
      "-$hasTshiftString-reflect.fits"$,
      "-$hasTshiftString-reflectPrim.fits"$,
      "-$noTshiftString-direct.fits"$
     ];
   
   print(fname);
   
   variable pl = xfig_plot_new(18,12);
   
   pl.world(0.1,50,1e-9,1;loglog);
   
   pl.title("\noindent Red: Irradiating Returning Spectrum; Yellow: Direct Spectrum \\"R+
+"Blue Xillver Reflection (blue) with its Primary Irrad (dashed) \\"R
+"Green: Xillver for Tin*1.4 (green) and its primary spectrum (dashed) "R);
   
   variable ii, n = length(fname);
   variable dat = Struct_Type[n];
 
   variable xillNorm = Double_Type[n];

   
   variable col_ind = [0,1,1,3,3,6];
   variable line_arr = [0,0,1,0,1,0];
   variable line_width=[3,1,1,1,1,3];
   
   _for ii(0,n-1){
      dat[ii] = fits_read_table(fname[ii]);
      
      
      variable en = 0.5*(dat[ii].bin_lo+dat[ii].bin_hi);
      pl.plot(en, dat[ii].flux;width=line_width[ii], color=CB_COLOR_SCHEME_NB[col_ind[ii]], line=line_arr[ii]);
      
      xillNorm[ii] = get_xillver_norm(en, dat[ii].flux);
      
      vmessage(" norm flux wrt xillver %s: %e", fname[ii], xillNorm[ii]);
   }

   pl.xlabel("Energy [keV]");
   pl.ylabel(" Counts / Bins ");

   variable rlo, rhi;
   (rlo,rhi) = getRgridFromHeader(fname[ii]);
   pl.add_object(xfig_new_text(sprintf("$T_\mathrm{in}=1$\,keV;  Radial Zone %.2f-%.2f$\,R_\mathrm{g}$"R,rlo,rhi)),
			       0.96,0.96,0.5,0.5;world0);
   

   variable retVal = EXIT_SUCCESS;
   
   if ( abs(xillNorm[0]-xillNorm[2]) > 1e-6){
      retVal = EXIT_SUCCESS;
      vmessage(" *** warning : incident and primary xillver spectrum do not have the same normaliztion [%e,%e]",
	      xillNorm[0],xillNorm[2]);
   }
   
   if ( xillNorm[1]>=xillNorm[2]){
      retVal = EXIT_SUCCESS;
      vmessage(" *** warning : reflected xillver norm is larger than primary xillver norm [%e,%e]",
	      xillNorm[1],xillNorm[2]);
   }      
   
   return (retVal,pl);
}

define testXillverPrimaryNormalization(){
   
   variable izones = [10,20,30,40];
   variable ii, n = length(izones);
   
   variable pl = Struct_Type[n];
   variable retVal = Int_Type[n];
   
   _for ii(0,n-1){
      (retVal[ii], pl[ii]) = testSingleXillverPrimaryNormalization(izones[ii]);
   }


   xfig_multiplot(pl;cols=2).render(dir_plots+"testXillverPrimaryNormalization.pdf";Verbose=-1);
   
   
   if (sum(retVal) == n*EXIT_SUCCESS){
      return EXIT_SUCCESS;
   } else {
      return EXIT_FAILURE;
   }

}


if (length(__argv)>0){
   () = testXillverPrimaryNormalization();
}
