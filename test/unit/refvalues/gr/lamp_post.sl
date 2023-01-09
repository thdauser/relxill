require("isisscripts");

define ut_h (h, a) {
  return sqrt((h^2 + a^2)/(h^2 - 2*h + a^2)); 
}



define gi(h,r,a)
{
   return ut_d(r,a)/ut_h(h,a);
}


define ut_h_v(h, a, v)
{
   variable sig = h^2 + a^2;
   variable del = sig - 2*h;

   return sqrt(  (del*sig) / (del^2 - v^2 * sig^2)  );
}



define lp_emisProfile_nonrelat(h,r){
   
   variable emis =  ( 1. / ( (r/h)^2 + 1 ) )^(3./2);
   emis /= (2*PI*h^2);

   return emis;
}


define lp_fluxBoost_obs2primary(h,a, Gamma){
   
      variable energShift = (kerr_lp_redshift(h,a)+1.);
   
      variable fluxBoost = energShift^Gamma;
   
      return fluxBoost;
   
}



%%% --------------------------- %%%
%%% DEPRECATED %%% 
require("isisscripts");
define gravRedshift(height, a){
   message(" gravRedshift DEPRECATED (now in the isisscripts!!), call with:  kerr_lp_redshift(height, a)");
   return kerr_lp_redshift(height, a);  
}
