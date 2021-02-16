require("isisscripts");
require("subs_returnrad.sl");

variable spin = 0.998;

define eval_kerrbb_model(lo,hi){
   fit_fun("kerrbb");
   set_par("kerrbb(1).rflag",qualifier("refl_frac",0));
   set_par(10,0);
   set_par(4,60);
   
   set_par(5,5);
   set_par(6,0.46);
   
   __set_hard_limits("kerrbb(1).hd",0.5,10);
   set_par(8,1./1.26,0,0.5,10);  %% T_eff = 1.26*fcol*m^-1/4*mdot^1/4

   set_par("kerrbb(1).a",spin);
   
   return eval_fun_keV(lo,hi)/(hi-lo) * 0.5*(lo+hi);
}

define get_2d_plot(fdat){

   variable pl = xfig_plot_new(13,10);

   pl.xlabel("Energy [keV]");
   pl.ylabel("Energy Flux $F_\nu$"R);
   pl.world(0.1,100,1e-6,7000;loglog);
   
   fdat = [fdat];
   variable m = length(fdat);
   
   variable dat = Struct_Type[m];

   variable ii, jj;
   _for jj(0,m-1){
      dat[jj] = get_2d_data(fdat[jj]);
   }


   
   _for jj(0,m-1){

      ccount=1; %% reset colour
      
      
      % make life easier
      variable elo = dat[jj].elo;
      variable ehi = dat[jj].ehi;
      variable emn = 0.5*(elo+ehi);
      variable efac = emn/(ehi-elo);
      
      
      pl.plot(emn,dat[jj].sumspec*efac; color="black",width=3,line=jj);
                  
      variable n = length(dat[jj].dat.spec[*,0]);
      variable mod_fac=9;
      _for ii(0,n-1){
	 
	 if (pl!=NULL){
	    if (ii mod mod_fac == 0) {
%	       pl.plot(emn,dat[jj].dat.spec[ii,*]*efac*dat[jj].area_ring[ii];color=get_col(), line=jj);
	       pl.plot(emn,dat[jj].dat.spec[ii,*]*efac;color=get_col(), line=jj);
	       
	       if(jj==0){
		  pl.add_object(xfig_new_text(sprintf("$r = %.2f\,r_\mathrm{g}$"R,
						      0.5*(dat[jj].dat.rlo[ii]+dat[jj].dat.rhi[ii]));
					      color=get_col(;same)),
				.78,0.95,-0.5,-0.5+1.1*ccount;world0);
	       }
	    }            	    
	 }
      }
   }

   pl.plot(emn,(dat[0].sumspec+dat[1].sumspec)*efac; color="red",width=3,depth=-10);


   if (qualifier_exists("plot_kerrbb")){
      variable val_kerrbb = eval_kerrbb_model(elo,ehi);
      variable kerrbb_fac = sum(dat[1].sumspec*efac*0.48) / sum(val_kerrbb);
      val_kerrbb *= kerrbb_fac;
      pl.plot(emn,val_kerrbb;width=3, color="blue");
      
      variable val_kerrbb_return = eval_kerrbb_model(elo,ehi;refl_frac=1);
      list_par;
      val_kerrbb_return *= kerrbb_fac;
      pl.plot(emn,val_kerrbb;width=3, color="blue", line=1);
   }
      
   return pl;
}


variable frel_dat = dir+[
			 "debug-testrr-spec2d-rr-bbody-reflect.fits",
%			 "debug-testrr-rframe-spec2d-rr-bbody.fits"
%			 "debug-testrr-rframe-spec2d-rr-bbody-xillver.fits",
			 "debug-testrr-spec2d-prim-bbody.fits",
			 "debug-testrr-rframe-spec2d-rr-bbody-xillver-prim.fits",
			 "debug-testrr-spec2d-rr-bbody.fits"
			];
variable pl2d = get_2d_plot(frel_dat);


pl2d.render("test_returnrad_relatSmearing.pdf");
