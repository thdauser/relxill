% -*- mode: slang; mode: fold -*-
require("isisscripts");
require("scripts/subs_filenames.sl");


%% xillver normalization N_xill of a spectrum (photons/bin), as defined for the xillver tables (Dauser+2016) and
%% calculated in the model (calc_xillver_norm_from_primary_spectrum): energy flux in 0.1-1000 keV relative to the
%% reference illumination of the table, 1e15/(4 pi)
variable global_norm_xillver = 1e15 / 4.0 / PI;

define get_xillver_norm(elo, ehi, spec){
   variable keV2erg = 1.602177e-09;
   variable ind = where(elo >= 0.1 and ehi <= 1000);
   return sum(spec[ind] * 0.5*(elo[ind]+ehi[ind])) * 1e20 * keV2erg / global_norm_xillver;
}


%% loads the spectra of the zones (photons/bin; rest frame: per unit area of the zone)
%% qualifier normfile: rest-frame returning spectra of the zones (rframe-bbodyRet); the spectrum of each zone is divided
%%    by N_xill of the returning spectrum of this zone (as in the model, where N_xill(specRet) scales the xillver table),
%%    i.e., all rest-frame spectra of a zone are normalized by the same factor,
%%    such that their relative normalization is kept and all zones are equally visible (see "Diagnostic output files"
%%    in doc/returnrad_disk_emission.md)
define get_2d_data(fname){
   vmessage(" Loading %s", fname);

   variable fd = fits_read_table(fname);
   variable elo = fd.ener[0,*];
   variable ehi = make_hi_grid(elo);
   variable ii, n = length(fd.rlo);

   variable normfile = qualifier("normfile", NULL);
   if (normfile != NULL){
      variable fnorm = fits_read_table(normfile);
      if (length(fnorm.rlo) != n){
         vmessage(" *** error: %s and %s have a different number of zones", fname, normfile);
         exit(1);
      }
      _for ii(0,n-1){
         fd.spec[ii,*] /= get_xillver_norm(elo, ehi, fnorm.spec[ii,*]);
      }
   }

   %% sum of the zones (only meaningful for the observed spectra)
   variable sumspec = elo*0;
   _for ii(0,n-1){
      sumspec += fd.spec[ii,*];
   }

   return struct{elo=elo, ehi=ehi, sumspec=sumspec, dat=fd};
}


define get_param_struct(spin, Tin, logxi, logn){
   return struct{spin=spin, Tin=Tin, logxi=logxi, logn=logn, incl=qualifier("incl",40.0)};
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

define eval_kerrbb_model(lo,hi,spin){
%{{{
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

define get_2d_plot_data(dat,rframe,spin){

   if (length(rframe)!=length(dat)){
      vmessage(" *** error : rframe array needs to be the same length as the number of files ");
      exit;
   }
   
   variable pl = xfig_plot_new(13,10);

   pl.xlabel("Energy [keV]");
   pl.ylabel(qualifier("ylabel", "Energy Flux $\nu\;F_\nu$"R));
   pl.world(0.1,250,1e-4,qualifier("ymax",7000);loglog);
   
   variable ii, jj;
   variable m = length(dat);
   _for jj(0,m-1){

      ccount=1; %% reset colour
      
      
      % make life easier
      variable elo = dat[jj].elo;
      variable ehi = dat[jj].ehi;
      variable emn = 0.5*(elo+ehi);
      variable efac = emn^2/(ehi-elo);
      

      if (not qualifier_exists("noSumSpec")){
	 pl.plot(emn,dat[jj].sumspec*efac; color="black",width=3,line=jj);
      }

      pl.add_object(xfig_new_text(sprintf("\fbox{$a = %.3f$}"R,spin);
					     size="\footnotesize"R),
				.03,0.99,-0.5,0.7;world0);
      
      variable n = length(dat[jj].dat.spec[*,0]);
      variable mod_fac=9;
      _for ii(0,n-1){
	 
	 if (pl!=NULL){
	    if (ii mod mod_fac == 0) {
	       pl.plot(emn,dat[jj].dat.spec[ii,*]*efac;color=get_col(), line=jj);
	       
	       if(jj==0){
		  pl.add_object(xfig_new_text(sprintf("$r = %.1f\,r_\mathrm{g}$"R,
						      0.5*(dat[jj].dat.rlo[ii]+dat[jj].dat.rhi[ii]));
					      color=get_col(;same),
					     size="\footnotesize"R),
				.78,0.99,-0.5,-0.3+0.9*(ccount-1);world0);
	       }
	    }            	    
	 }
      }
      if (not qualifier_exists("noSumSpec")){
	 vmessage(" summed energy flux: %.4e (frac: %.3f)", sum(dat[jj].sumspec*elo),
		  sum(dat[jj].sumspec*elo)/sum(dat[0].sumspec*elo));
      }
   }

   if (not qualifier_exists("noSumSpec")){
      pl.plot(emn,(dat[0].sumspec+dat[1].sumspec)*efac; color="red",width=3,depth=-10);
   }
  
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

%% rframe[jj]=1: rest-frame file, normalized per zone by N_xill of the rest-frame returning spectrum of the zone (sums over
%% the zones are not meaningful then, use the qualifier noSumSpec)
define get_2d_plot(fdat,rframe,spin){

   fdat = [fdat];
   variable m = length(fdat);
   
   variable dat = Struct_Type[m];
   variable normfile = getFilenameStruct().rframe.specRet;

   variable jj;
   _for jj(0,m-1){
      if (rframe[jj]==1){
	 dat[jj] = get_2d_data(fdat[jj]; normfile=normfile);
      } else {
	 dat[jj] = get_2d_data(fdat[jj]);
      }
   }

   variable qual = __qualifiers();
   if (any(rframe==1) && not qualifier_exists("ylabel")){
      qual = struct_combine(qual, struct{ylabel="$\nu\;F_\nu$ / $N_\mathrm{xill}$(returning spectrum of the zone)"R});
   }

   return get_2d_plot_data(dat, rframe, spin;; qual);
}

