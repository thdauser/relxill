require("isisscripts");
require("subs_returnrad.sl");

variable spin = 0.0;

define eval_kerrbb_model(lo,hi){
   fit_fun("kerrbb");
   set_par(9,0);
   set_par(10,0);
   set_par(4,60);
   __set_hard_limits("kerrbb(1).hd",0.5,10);
   set_par(8,1./1.26,0,0.5,10);  %% T_eff = 1.26*fcol*m^-1/4*mdot^1/4

   set_par("kerrbb(1).a",spin);
   
   return eval_fun_keV(lo,hi)/(hi-lo) * 0.5*(lo+hi);
}

define get_2d_plot(fdat_rframe,fdat_rel){

   variable pl = xfig_plot_new(13,10);
   variable plf = xfig_plot_new(13,6);

   pl.xlabel("Energy [keV]");
   pl.ylabel("Energy Flux $F_\nu$"R);
   pl.world(0.00035,1000,1e-8,1000;loglog);

   plf.world(1,1000,0.1,1000;loglog);
   
   variable m = 2;
   variable dat = Struct_Type[m];
   dat[0] = get_2d_data(fdat_rframe;diskarea);
   dat[1] = get_2d_data(fdat_rel);

   variable flux_single_zones = Array_Type[m];
   
   
   
   variable ii, jj;
   _for jj(0,m-1){

      ccount=1; %% reset colour
      
      
      % make life easier
      variable elo = dat[jj].elo;
      variable ehi = dat[jj].ehi;
      variable emn = 0.5*(elo+ehi);
      variable efac = emn/(ehi-elo);

      efac=1.;
      
      pl.plot(emn,dat[jj].sumspec*efac; color="black",width=3,line=jj);
                  
      variable n = length(dat[jj].dat.spec[*,0]);

      flux_single_zones[jj]  = Double_Type[n];

      variable mod_fac=8;
      _for ii(0,n-1){

	 flux_single_zones[jj][ii] = sum(dat[jj].dat.spec[ii,where(0.01<=elo<1000.)]);
	 
	 if (pl!=NULL){
	    if (ii mod mod_fac == 0) {
	       pl.plot(emn,dat[jj].dat.spec[ii,*]*efac;color=get_col(), line=jj);
	       print(sum(dat[jj].dat.spec[ii,*]));
	       
	       if(jj==1){
		  pl.add_object(xfig_new_text(sprintf("$r = %.2f\,r_\mathrm{g}$"R,
						      0.5*(dat[jj].dat.rlo[ii]+dat[jj].dat.rhi[ii]));
					      color=get_col(;same)),
				.78,0.95,-0.5,-0.5+1.1*ccount;world0);
	       }
	    }            	    
	 }
      }
   }


   % pl.plot(emn,(dat[0].sumspec+dat[1].sumspec)*efac; color="red",width=3,depth=-10);


   
   variable val_kerrbb = eval_kerrbb_model(elo,ehi);
   val_kerrbb *= sum(dat[1].sumspec*efac) / sum(val_kerrbb);
%   pl.plot(emn,val_kerrbb;width=3, color="blue");

   
   
   _for jj(0,m-1){
      plf.plot(0.5*(dat[jj].dat.rlo+dat[jj].dat.rhi),
	       flux_single_zones[jj];
	       line=jj,width=3);
   }
   plf.plot(0.5*(dat[0].dat.rlo+dat[0].dat.rhi),
	    flux_single_zones[0]/flux_single_zones[1];width=3,color=1);
   
   
   return xfig_new_vbox_compound(pl,plf);
}


%variable frel_dat = dir+[spec2d-rr-bbody.fits","debug-testrr-spec2d-prim-bbody.fits"];
%variable frel_dat = dir+["debug-testrr-spec2d-rr-bbody.fits","testrr-rframe-rr-bbody.fits"];
variable frel_dat = dir+["debug-testrr-rframe-spec2d-prim-bbody.fits","debug-testrr-spec2d-prim-bbody.fits"];

variable pl2d = get_2d_plot(frel_dat[0],frel_dat[1]);


pl2d.render("test_returnrad_primaryRad.pdf");
