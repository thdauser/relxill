% -*- mode: slang; mode: fold -*-
require("isisscripts");
require("scripts/subs_basic_equations.sl");


define get_2d_data(fname, spin){ %{{{
   
   vmessage(" Loading %s (assuming spin=%.3f)", fname, spin);
      
   variable fd = fits_read_table(fname);
   
   variable ii;
   variable n = length(fd.rlo);
   variable nener = length(fd.spec[0,*]);
   
   variable elo = fd.ener[0,*];
   variable ehi = make_hi_grid(elo);

   variable emn = 0.5*(elo+ehi);

   variable area_ring = calc_proper_area_ring_relxill(fd.rlo,fd.rhi,spin);

   variable sumspec = elo*0;
   _for ii(0,n-1){

      if (qualifier("diskarea",0)==1){
	 fd.spec[ii,*] *= area_ring[ii];
      }
      sumspec += fd.spec[ii,*];
      
   }   
   
   return struct{elo=elo, ehi=ehi, sumspec=sumspec, area_ring=area_ring, dat=fd};
}
%}}}



define get_xillver_norm(en,flux){
   

   variable ind = where(0.1 < en < 1000);
   
   variable norm_xill = 1e15 / 4.0 / PI;
   variable keV2erg = 1.602177e-09;
      
   variable reflux = en * 1e20 * keV2erg * flux;
   
   variable norm_flux = sum(reflux[ind]);
   
   return norm_flux / norm_xill;
}


variable ccount=0;
define get_col(){
   if (qualifier_exists("same")) {
      ccount--;
   }
      variable col = [CB_COLOR_SCHEME_NB,CB_COLOR_SCHEME2_NB][ccount mod 14];
   ccount++;
   return col;
}


define eval_kerrbb_model(lo,hi,spin){ %{{{
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
%}}}


define get_2d_plot(fdat,rframe,spin){

   if (length(rframe)!=length(fdat)){
      vmessage(" *** error : rframe array needs to be the same length as the number of files ");
      exit;
   }
   
   variable pl = xfig_plot_new(13,10);

   pl.xlabel("Energy [keV]");
   pl.ylabel("Energy Flux $F_\nu$"R);
   pl.world(0.1,30,1e-4,7000;loglog);
   
   fdat = [fdat];
   variable m = length(fdat);
   
   variable dat = Struct_Type[m];

   variable ii, jj;
   _for jj(0,m-1){
      dat[jj] = get_2d_data(fdat[jj],spin;diskarea=rframe[jj]);
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

   if (not qualifier_exists("noSumSpec")){
      pl.plot(emn,(dat[0].sumspec+dat[1].sumspec)*efac; color="red",width=3,depth=-10);
   }

   print(sum(dat[0].sumspec*elo));
   print(sum(dat[1].sumspec*elo));
   print(sum(dat[2].sumspec*elo));
  
   if (qualifier_exists("title")){
      pl.title(qualifier("title","description missing");size="\footnotesize"R);
   }

   if (qualifier_exists("plot_kerrbb")){
      variable val_kerrbb = eval_kerrbb_model(elo,ehi,spin);
      variable ind = where (elo>0.1);
      variable kerrbb_fac = sum(dat[1].sumspec[ind]*efac[ind]) / sum(val_kerrbb[ind]);
      val_kerrbb *= kerrbb_fac;
      pl.plot(emn,val_kerrbb;width=3, color="blue");
      
      variable val_kerrbb_return = eval_kerrbb_model(elo,ehi, spin;refl_frac=1);
      list_par;
      val_kerrbb_return *= kerrbb_fac;
      pl.plot(emn,val_kerrbb;width=3, color="blue", line=1);
   }
      
   return pl;
}
