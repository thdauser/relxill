require("isisscripts");
require("scripts/subs_testSetup.sl");
require("scripts/subs_returnrad.sl");


define testRframeReturnradRefdata(){

   
   variable spin = 0.998;
   
   variable dir = "build/";
   variable fname_p = dir+"testrr-rframe-prim-bbody.fits";
   variable fname_r = dir+"testrr-rframe-rr-bbody.fits";
   %variable fname_p = dir+"debug-testrr-bbody-rframe-specPri.fits";
   %variable fname_r = dir+"debug-testrr-bbody-rframe-specRet.fits";
   
   variable dp = get_2d_data(fname_p,spin;diskarea=1);
   variable dr = get_2d_data(fname_r,spin;diskarea=1);
   
   % make life easier
   variable elo = dp.elo;
   variable ehi = dp.ehi;
   variable emn = 0.5*(elo+ehi);
   
   
   variable pl = xfig_plot_new(17,10);
   variable plr = xfig_plot_new(17,5);
   
   %pl.world(0.1,20,1e-1,5000;loglog);
   
   %pl.world(0.1,60,2e-2,340;loglog);
   pl.world(0.01,50,1e-1,400;loglog);
   plr.world(0.01,50,0.99,1.02;xlog);
   
   pl.xlabel("Energy [keV]");
   pl.ylabel("Energy Flux $F_\nu$"R);
   plr.xlabel("Energy [keV]");
   plr.ylabel("Ratio"R);
   
     
   variable fref = dir_refdata+"refvalues-rrbbody.fits";
   variable ref = fits_read_table(fref);

   
   variable testTitle = "\noindent Comparing Rest-Frame Prim (dashed) and Ret (solid) with Refdata"R;
   testTitle += sprintf(" %s ",fref );

   
   variable efac = emn/(ehi-elo);
   %%% difference in normalization (calculating the overall disk area)
   %%% between relline and the return-rad internal calculation 
   efac *= 4;  
   dp.sumspec *= efac;
   dr.sumspec *= efac;
   
   variable pcol = "blue";
   variable rcol = "red";
   
   pl.plot(emn,dp.sumspec; color=pcol,width=3,line=1);
   
   pl.plot(emn,dr.sumspec; color=rcol,width=3,sym="point");
   
   pl.plot(emn,dr.sumspec+dp.sumspec; color="black",width=3);
   
   
   variable ii, n = length(dp.dat.spec[*,0]);
   variable mod_fac=8;
   _for ii(0,n-1){
      
      if (pl!=NULL){
	 if (ii mod mod_fac == 0) {
	    pl.plot(emn,dr.dat.spec[ii,*]*efac;color=get_col());
	    pl.plot(emn,dp.dat.spec[ii,*]*efac;color=get_col(;same),line=1);
	    pl.add_object(xfig_new_text(sprintf("$r = %.2f\,r_\mathrm{g}$"R,
						0.5*(dp.dat.rlo[ii]+dp.dat.rhi[ii]));
					color=get_col(;same)),
			  .78,0.95,-0.5,-0.5+1.1*ccount;world0);
	 }            	    
      }
   }
   
   
   variable eref = 0.5*(ref.elo+ref.ehi);
   ref.spec_ret *= eref/(ref.ehi-ref.elo);
   ref.spec_pri *= eref/(ref.ehi-ref.elo);
   
   %% %REF Plot
   pl.plot(eref,ref.spec_ret;color="green3",line=1);
   pl.plot(eref,ref.spec_pri;color="green3",line=1);
   
   
   variable rebin_r_sumspec = interpol_points(eref,emn,dr.sumspec);
   variable rebin_p_sumspec = interpol_points(eref,emn,dp.sumspec);
   
   plr.plot(eref,rebin_r_sumspec/ref.spec_ret;color=rcol,width=3);
   plr.plot(eref,rebin_p_sumspec/ref.spec_pri;color=pcol,line=1,width=3);
   
   %plr.plot(emn,dr.sumspec/(ref.spec_ret);color=rcol,width=3);
   %plr.plot(emn,dp.sumspec/(ref.spec_pri);color=pcol,line=1,width=3);
   
   
   pl.title(testTitle);
   
   
   mkdir(dir_plots);
   xfig_multiplot(pl,plr).render(dir_plots+"testRframeReturnradRefvalues.pdf");
   

   variable prec = 1e-3;
   variable isPrimaryApproxEq = ( abs(mean(rebin_p_sumspec/ref.spec_pri)-1) < prec )  ? 1 : 0;
   variable isReturnApproxEq = ( abs(mean(rebin_r_sumspec/ref.spec_ret)-1) < prec )  ? 1 : 0;

   if (isPrimaryApproxEq==1 && isReturnApproxEq==1) {
      return EXIT_SUCCESS;
   } else {
      vmessage(" *** error : restframe spectra do not agre with reference data ");
      return EXIT_FAILURE;
   }
   
}
