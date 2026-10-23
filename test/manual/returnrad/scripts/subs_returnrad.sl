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

%% plots the spectra of every 9th zone of the data dat (one line style per file)
%% qualifiers:
%%   lines      : line style of each file (0: solid, 1: dashed, 2: dotted; default: the index of the file)
%%   noSumSpec  : do not plot the sum over the zones of each file (black, same line style)
%%   total      : indices of the files whose sums over the zones are added up and plotted in red (e.g., direct +
%%                reflected = model)
%%   ylabel, ymax, title
define get_2d_plot_data(dat,spin){

   variable pl = xfig_plot_new(13,10);

   pl.xlabel("Energy [keV]");
   pl.ylabel(qualifier("ylabel", "Energy Flux $\nu\;F_\nu$"R));
   pl.world(0.1,250,1e-4,qualifier("ymax",7000);loglog);

   pl.add_object(xfig_new_text(sprintf("\fbox{$a = %.3f$}"R,spin); size="\footnotesize"R),
		 .03,0.99,-0.5,0.7;world0);

   variable ii, jj;
   variable m = length(dat);
   variable lines = qualifier("lines", [0:m-1]);
   variable emn, efac;
   _for jj(0,m-1){

      ccount=1; %% reset colour

      variable elo = dat[jj].elo;
      variable ehi = dat[jj].ehi;
      emn = 0.5*(elo+ehi);
      efac = emn^2/(ehi-elo);

      if (not qualifier_exists("noSumSpec")){
	 pl.plot(emn,dat[jj].sumspec*efac; color="black",width=3,line=lines[jj]);
	 vmessage(" summed energy flux of file %d: %.4e (relative to file 0: %.4f)", jj, sum(dat[jj].sumspec*emn),
		  sum(dat[jj].sumspec*emn)/sum(dat[0].sumspec*emn));
      }

      variable n = length(dat[jj].dat.spec[*,0]);
      variable mod_fac=9;
      _for ii(0,n-1){
	 if (ii mod mod_fac == 0) {
	    pl.plot(emn,dat[jj].dat.spec[ii,*]*efac;color=get_col(), line=lines[jj]);

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

   variable itot = qualifier("total", NULL);
   if (itot != NULL){
      variable sum_total = dat[itot[0]].sumspec*0;
      foreach jj (itot){
	 sum_total += dat[jj].sumspec;
      }
      pl.plot(emn, sum_total*efac; color="red",width=3,depth=-10);
   }

   if (qualifier_exists("title")){
      pl.title(qualifier("title","description missing");size="\footnotesize"R);
   }

   return pl;
}

%% plot of the spectra of the zones of the files fdat (see get_2d_plot_data for the qualifiers)
%% qualifier rframe: rest-frame files, normalized per zone by N_xill of the rest-frame returning spectrum of the zone
%%    (the sums over the zones are not meaningful then and not plotted)
define get_2d_plot(fdat,spin){

   fdat = [fdat];
   variable m = length(fdat);
   variable dat = Struct_Type[m];
   variable qual = __qualifiers();

   variable jj;
   _for jj(0,m-1){
      if (qualifier_exists("rframe")){
	 dat[jj] = get_2d_data(fdat[jj]; normfile=getFilenameStruct().rframe.specRet);
      } else {
	 dat[jj] = get_2d_data(fdat[jj]);
      }
   }

   if (qualifier_exists("rframe")){
      qual = struct_combine(qual, struct{noSumSpec});
      if (not qualifier_exists("ylabel")){
	 qual = struct_combine(qual, struct{ylabel="$\nu\;F_\nu$ / $N_\mathrm{xill}$(returning spectrum of the zone)"R});
      }
   }

   return get_2d_plot_data(dat, spin;; qual);
}

