require("isisscripts");

  

variable x0_lo = [1,2,4,5,6,7,9];
variable x0_hi = [2,4,5,6,7,9,10];

variable y0 = [1,4,3,2,4,1,2];

variable xn_lo = [1,3,4,8.5];
variable xn_hi = [3,4,8.5,10];

variable x0_m = 0.5*(x0_lo+x0_hi);
variable xn_m = 0.5*(xn_lo+xn_hi);

print(interpol_points(xn_m,x0_m,y0));

print(rebin_mean(xn_lo,xn_hi,x0_lo,x0_hi,y0));

