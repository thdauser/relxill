require("isisscripts");

evalfile("load_model.sl");

variable ff = "xillver";
variable param = "z";
%variable param_vals = [0,0.001,0.01,0.1,0.2,0.5];
variable param_vals = [1.0];


%%%%%%%%%%

variable lo, hi;
(lo,hi) = log_grid(0.1,800,500);

fit_fun(ff);
xlin;ylog;

xrange(1.1,1.2);

variable dat = fits_read_table("refdata_localModels/xillver/xillver_random_refdat_0030.fits_comparison.tmp");
hplot(dat.bin_lo,dat.bin_hi, dat.model / (dat.bin_hi - dat.bin_lo ));
ohplot(dat.bin_lo,dat.bin_hi, dat.value/ (dat.bin_hi - dat.bin_lo ));

load_par("refdata_localModels/xillver/xillver_random_refdat_0030.fits.par");
%variable val = eval_fun_keV(dat.bin_lo,dat.bin_hi);
%ohplot(dat.bin_lo,dat.bin_hi,val);

variable zpar = get_par("*.z");

variable ii;
_for ii(0,length(param_vals)-1) {
   set_par("*."+param,zpar*param_vals[ii] );
   
%   variable val = eval_fun_keV(dat.bin_lo,dat.bin_hi);
%   ohplot(dat.bin_lo,dat.bin_hi,val);
   variable val = eval_fun_keV(lo,hi)/(hi-lo);
   ohplot(lo,hi,val);
   
}


%variable val = eval_fun_keV(lo,hi);




