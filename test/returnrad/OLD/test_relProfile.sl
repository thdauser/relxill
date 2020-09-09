require("isisscripts");
require("subs_returnrad.sl");

%variable dat = fits_read_table(dir+"debug-testrr-relline-profile.fits");
variable rprof = get_2d_data(dir+"debug-testrr-relline-profile.fits");

variable ii;
variable n = length(rprof.area_ring);
variable sumrad = Double_Type[n];
_for ii(0,n-1){
   sumrad[ii] = sum(rprof.dat.spec[ii,*]);
}

yrange(1e-4,1);
xrange(1,1e3);
xlog;ylog;
variable r = 0.5*(rprof.dat.rlo+rprof.dat.rhi);
plot(r, rprof.area_ring*r^(-3));
oplot(r, sumrad*r^(-3));

%ylin;
%yrange(0.7,1.3);
%plot(r, rprof.area_ring/sumrad);
