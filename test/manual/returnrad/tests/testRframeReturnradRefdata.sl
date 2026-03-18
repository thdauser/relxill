require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");


define model_eval(){

   putenv("RELXILL_WRITE_FILES=1");
   load_xspec_local_models("./build/");
   variable ff = "relxillBBxill";
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


define plot_single_ref_data(dp, dp_ref){
%%%%%
   variable pl = xfig_plot_new(17,10);
   
   pl.world(0.01,150,1e-1,400;loglog);
   
   pl.xlabel("Energy [keV]");
   pl.ylabel("Energy Flux $F_\nu$"R);


   % make life easier
   variable emn = 0.5*(dp.elo+dp.ehi);
   variable emn_ref = 0.5*(dp_ref.elo+dp_ref.ehi);

   variable efac = emn/(dp.ehi-dp.elo);
   variable efac_ref = emn_ref/(dp_ref.ehi-dp_ref.elo);

   variable pcol = "blue";
   variable rcol = "red";
   pl.plot(emn_ref,dp_ref.sumspec*efac_ref; color="black",width=3,line=1);
   pl.plot(emn,dp.sumspec*efac_ref; color="red",width=3);

   % xlog;ylog;
   % yrange(0.1,400);
   % plot(emn,dp.sumspec*efac);
   % sleep(10);
   
%   pl.plot(emn,dr.sumspec; color=rcol,width=3,sym="point");   
%   pl.plot(emn,dr.sumspec+dp.sumspec; color="black",width=3);

   variable ii, n = length(dp.dat.spec[*,0]);
   variable mod_fac=8;
   ccount = 0;
   _for ii(0,n-1){
      
      if (pl!=NULL){
	 if (ii mod mod_fac == 0) {
	    pl.plot(emn,dp_ref.dat.spec[ii,*]*efac_ref;color="black",line=1);
	    pl.plot(emn,dp.dat.spec[ii,*]*efac;color=get_col(),line=0);
	    pl.add_object(xfig_new_text(sprintf("$r = %.2f\,r_\mathrm{g}$"R,
						0.5*(dp.dat.rlo[ii]+dp.dat.rhi[ii]));
					color=get_col(;same)),
			  .78,0.95,-0.5,-0.5+1.1*ccount;world0);
	 }            	    
      }
   }

   return pl;
}


define testRframeReturnradRefdata(){
%%%%%
   variable spin = model_eval();

   variable fnames = get_fnames();

   variable dp = get_2d_data(fnames.prim,spin;diskarea=1);
   variable dr = get_2d_data(fnames.refl,spin;diskarea=1);

   variable dir_ref = "refdata/";
   variable dp_ref = get_2d_data(dir_ref+fnames.prim,spin;diskarea=1);
   variable dr_ref = get_2d_data(dir_ref+fnames.refl,spin;diskarea=1);


   
   variable prec = 1e-5;
   variable isPrimaryApproxEq = is_approx_equal(dp.sumspec, dp_ref.sumspec, prec );
   variable isReturnApproxEq = is_approx_equal(dr.sumspec, dr_ref.sumspec, prec );

   variable status = EXIT_SUCCESS;
   if (isPrimaryApproxEq==0 && isReturnApproxEq==0) {
      vmessage(" *** error : restframe spectra do not agre with reference data ");
      return EXIT_FAILURE;
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

