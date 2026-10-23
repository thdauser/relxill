require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");


define model_eval(){

   putenv("RELXILL_WRITE_FILES=1");
   load_xspec_local_models("./build/");
   variable ff = "relxillBB";
   variable spin = 0.998;

   fit_fun(ff);
   set_par("*.a",spin);

   () = eval_fun(1,2);

   return spin;
}

define get_fnames(){

   variable f = getFilenameStruct();

   return struct{prim=f.rframe.specPri , refl=f.rframe.specRet};
}


%% rest-frame spectra of every 8th zone (solid) and of the reference data (dashed), normalized per zone by N_xill of
%% the reference returning spectrum
define plot_single_ref_data(dp, dp_ref){
%%%%%
   variable pl = xfig_plot_new(17,10);
   
   pl.world(0.01,250,1e-4,7000;loglog);
   
   pl.xlabel("Energy [keV]");
   pl.ylabel("$\nu\;F_\nu$ / $N_\mathrm{xill}$(reference returning spectrum of the zone)"R);

   variable emn = 0.5*(dp.elo+dp.ehi);
   variable efac = emn^2/(dp.ehi-dp.elo);
   variable emn_ref = 0.5*(dp_ref.elo+dp_ref.ehi);
   variable efac_ref = emn_ref^2/(dp_ref.ehi-dp_ref.elo);

   variable ii, n = length(dp.dat.spec[*,0]);
   variable mod_fac=8;
   ccount = 0;
   _for ii(0,n-1){
      if (ii mod mod_fac == 0) {
	 pl.plot(emn_ref,dp_ref.dat.spec[ii,*]*efac_ref;color="black",line=1);
	 pl.plot(emn,dp.dat.spec[ii,*]*efac;color=get_col(),line=0);
	 pl.add_object(xfig_new_text(sprintf("$r = %.2f\,r_\mathrm{g}$"R,
					     0.5*(dp.dat.rlo[ii]+dp.dat.rhi[ii]));
				     color=get_col(;same)),
		       .78,0.95,-0.5,-0.5+1.1*ccount;world0);
      }
   }

   return pl;
}


%% every zone agrees with the reference data
define is_approx_equal_zones(d, d_ref, prec){
   variable ii, n = length(d_ref.dat.rlo);
   if (length(d.dat.rlo) != n) return 0;
   _for ii(0,n-1){
      if (is_approx_equal(d.dat.spec[ii,*], d_ref.dat.spec[ii,*], prec) == 0){
	 vmessage("   zone %d (r = %.2f-%.2f) deviates from the reference data", ii, d.dat.rlo[ii], d.dat.rhi[ii]);
	 return 0;
      }
   }
   return 1;
}


define testRframeReturnradRefdata(){
%%%%%
   variable spin = model_eval();

   variable fnames = get_fnames();

   %% the rest-frame spectra of the zones, all normalized per zone by N_xill of the reference returning spectrum of the zone (the
   %% same factor for the test and the reference data, so any change of the normalization is detected)
   variable dir_ref = "refdata/";
   variable normfile = dir_ref+fnames.refl;
   variable dp = get_2d_data(fnames.prim; normfile=normfile);
   variable dr = get_2d_data(fnames.refl; normfile=normfile);

   variable dp_ref = get_2d_data(dir_ref+fnames.prim; normfile=normfile);
   variable dr_ref = get_2d_data(dir_ref+fnames.refl; normfile=normfile);

   variable prec = 1e-3;
   variable isPrimaryApproxEq = is_approx_equal_zones(dp, dp_ref, prec);
   variable isReturnApproxEq = is_approx_equal_zones(dr, dr_ref, prec);

   variable status = EXIT_SUCCESS;
   if (isPrimaryApproxEq==0 || isReturnApproxEq==0) {
      vmessage(" *** error : restframe spectra do notwith reference data \n\n");
      status = EXIT_FAILURE;
   } else {
      vmessage(" successful : restframe spectra agree with reference data \n\n");
   }


   %%% PLOT: Only in case of an error, or if qualifier ist given %%% 
   if (status==EXIT_FAILURE || qualifier_exists("plot")){
      
      variable testTitle = "\noindent Comparing Rest-Frame Prim (dashed, black) and Ret (solid, color) with Refdata"R;      
      
      variable pl_pri = plot_single_ref_data(dp, dp_ref);   
      variable pl_ret = plot_single_ref_data(dr, dr_ref);   
      
      mkdir(dir_plots);
      pl_pri.title(testTitle);
      xfig_multiplot(pl_pri, pl_ret).render(dir_plots+"testRframeReturnradRefvalues.pdf");
   }

   return status;
}

