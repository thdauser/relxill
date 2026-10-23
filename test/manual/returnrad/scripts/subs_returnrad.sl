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

%% plots the spectra of the zones of the data dat (one line style per data set)
%% qualifiers:
%%   nzones     : number of zones shown, evenly spaced in the zone index from the innermost zone to the zone with the
%%                mean radius closest to rmax (default 4)
%%   rmax       : radius of the outermost zone shown (default 100 Rg; the outer zones are cool and emit at very low
%%                energies)
%%   lines      : line style of each file (0: solid, 1: dashed, 2: dotted; default: the index of the file)
%%   noSum      : do not plot the sum over the zones of each data set (black, in its line style)
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

      if (not qualifier_exists("noSum")){
	 pl.plot(emn,dat[jj].sumspec*efac; color="black",width=3,line=lines[jj]);
	 vmessage(" energy flux integrated over the zones of data set %d: %.4e (relative to data set 0: %.4f)", jj,
		  sum(dat[jj].sumspec*emn), sum(dat[jj].sumspec*emn)/sum(dat[0].sumspec*emn));
      }

      variable n = length(dat[jj].dat.spec[*,0]);
      variable nzones = qualifier("nzones", 4);
      variable rmid = 0.5*(dat[jj].dat.rlo+dat[jj].dat.rhi);
      variable dist_rmax = abs(rmid - qualifier("rmax", 100.0));
      variable imax = where(dist_rmax == min(dist_rmax))[0];
      variable izones = nint([0:nzones-1] * imax / (nzones-1.0));
      foreach ii (izones){
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


   if (qualifier_exists("title")){
      pl.title(qualifier("title","description missing");size="\footnotesize"R);
   }

   return pl;
}

%% sum of the spectra of the zones of several files (same zones and energy grid)
define add_2d_data(d_a, d_b){
   variable fd = @d_a.dat;
   fd.spec = d_a.dat.spec + d_b.dat.spec;
   return struct{elo=d_a.elo, ehi=d_a.ehi, sumspec=d_a.sumspec + d_b.sumspec, dat=fd};
}

%% plot of the spectra of the zones (see get_2d_plot_data for the qualifiers); each element of fdat is one data set,
%% given by a file name or by an array of file names whose spectra are added up (e.g., {fdirect, [fdirect, frefl]})
%% qualifier rframe: rest-frame files, normalized per zone by N_xill of the rest-frame returning spectrum of the zone
%%    (the sums over the zones are not meaningful then and not plotted)
define get_2d_plot(fdat,spin){

   variable m = length(fdat);
   variable dat = Struct_Type[m];
   variable qual = __qualifiers();

   variable jj, fname;
   _for jj(0,m-1){
      dat[jj] = NULL;
      foreach fname ([fdat[jj]]){
	 variable d = qualifier_exists("rframe") ?
	    get_2d_data(fname; normfile=getFilenameStruct().rframe.specRet) : get_2d_data(fname);
	 dat[jj] = (dat[jj] == NULL) ? d : add_2d_data(dat[jj], d);
      }
   }

   if (qualifier_exists("rframe")){
      qual = struct_combine(qual, struct{noSum});
   }
   if (qualifier_exists("rframe") && not qualifier_exists("ylabel")){
      qual = struct_combine(qual, struct{ylabel="$\nu\;F_\nu$ / $N_\mathrm{xill}$(returning spectrum of the zone)"R});
   }

   return get_2d_plot_data(dat, spin;; qual);
}

