#!/usr/bin/env isis-script
% -*- mode: slang; mode: fold -*-

%% we can load any model
variable ff = "*";
if (length(__argv)>1){
   ff = __argv[1];
}
load_xspec_local_models("build/librelxill.so");

__set_hard_limits("relxilllp","h",-100,1000);


require("isisscripts");
require("fits_model_struct");


variable DATA_DIR = "refdata/";

variable sum_name = "relxill_sum.pdf";
variable glob_str = sprintf("test*_%s_*.fits",ff);
variable fnames = glob(DATA_DIR+sprintf(glob_str));
fnames = fnames[array_sort(fnames)];

define add_labels(pl,goodness,fn){
   
   pl.add_object(xfig_new_text(sprintf(`%s [goodness %.3e] --- %s`,
				       get_fit_fun(),goodness,
				       strreplace(fn,`_`,`\_`))),
		 0.5,1.05,0,0;world0);
   
   variable params = get_params();
   
   variable ii, n=length(params);
   
   variable str;
   _for ii(0,n-1,1){
      str = sprintf(`%s : %.3e`, strreplace(params[ii].name,`_`,`\_`),
		    params[ii].value);
      pl.add_object(xfig_new_text(str;
				 size=`\footnotesize`),
		    1.02,0.95,-0.5,0+0.8*ii;world0);            
   }
}

define check_single_model(fn){
   
   variable dat = fits_read_model_struct(fn);
   variable m_dat = eval_fun_keV(dat.bin_lo,dat.bin_hi);

   save_par(fn+".par");
   
   variable pl = xfig_plot_new(15,8);
   variable plr = xfig_plot_new(15,4);

   
   variable elo = dat.bin_lo;
   variable ehi = dat.bin_hi;
   
   variable ind_no0 = where(dat.value != 0);
   dat.value = dat.value[ind_no0];
   m_dat = m_dat[ind_no0];
   
   pl.hplot([elo,ehi[-1]], dat.value/(ehi-elo)*(0.5*(elo+ehi))^2;
	    loglog);
   pl.hplot([elo,ehi[-1]], m_dat/(ehi-elo)*(0.5*(elo+ehi))^2 ;
	    loglog, color="red", depth=-1 );
   
   plr.hplot([elo,ehi[-1]], dat.value*0+1; loglog);
   
   plr.hplot([elo,ehi[-1]], m_dat/dat.value; loglog, color="red");

   pl.ylabel(" Energy Flux ");
   plr.ylabel(" Ratio New/Old ");
   plr.xlabel("Energy [keV]" );

   
   variable ind = where(elo>0.2);
   variable goodness =  sqrt(sum((m_dat[ind]/dat.value[ind]-1)^2))/length(m_dat[ind]);
   
   vmessage(" GOODNESS (sqr-distance / num_bins ): %.2e ", goodness);

      add_labels(pl,goodness,fn);

   
   xfig_multiplot(pl,plr).render(fn+".pdf");
   
   return goodness;
}

variable fn;

variable ii;
variable n = length(fnames);

variable goodness = Double_Type[n];

_for ii(0,n-1,1){
   goodness[ii] = check_single_model(fnames[ii]);   
}

variable out = struct{name=fnames, goodness=goodness};
struct_filter(out, array_sort(goodness;dir=-1));


message(" \n **** SUMMARY **** \n");
_for ii(0,n-1,1){
   vmessage(" %s  \t %.3e ",out.name[ii], out.goodness[ii]);
}


() = system(sprintf("pdftk %s%s.pdf cat output %s%s",DATA_DIR,glob_str,DATA_DIR,sum_name));


#iffalse
variable lo, hi;
(lo, hi) = log_grid(0.9,10,20);

variable hist = histogram(goodness,lo,hi);
hplot(lo,hi,hist);