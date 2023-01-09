define q2(del, h, a) {
  return (sin(del)^2*(h^2+a^2)^2)/(h^2+a^2-2*h)-a^2;
}

define v_eta (a,q2,eta) {
  return -a^2*eta^4 + (a^2 - q2)*eta^2 + q2;
}


define v_r (a,q2,r){
  return r^4 + (a^2 - q2)*r^2 + 2*(a^2+q2)*r - a^2*q2;
}

define v_phi(a,r){
   
   variable vphi = r^2-2*a*sqrt(r)+a^2;
 
   variable delta = r^2-2*r+a^2;
   
   vphi /=  sqrt(delta) * (r^(3./2) + a) ;
   
   return vphi;
}


define ut_d (r, a) {
  return (
	  (r*sqrt(r)+a)/
	  (sqrt(r)*sqrt(r^2 -3*r + 2*a*sqrt(r)))
	 );
}

define ut_d_lambda (r, a, lambda) {
  return ((r*sqrt(r)+a - lambda )/(sqrt(r)*sqrt(r^2 -3*r + 2*a*sqrt(r))));
}

define del_inc(q,r,a)
{
   return q/(ut_d(r,a)*r);
}


define cal_gamma(r,a){
   variable v = (r^2-2*a*sqrt(r)+a^2)
     / (sqrt(r^2 - 2*r + a^2) * (r^(3./2.) + a) ) ;

   variable ga = 1./ sqrt(1.-v^2);
	
   return ga;
}

define cal_a_gr(r,a){
   variable rho2 = r^2 + a^2;  %% this is rho^2
   variable del = r^2 - 2*r + a^2;

   variable area_gr = 2*PI / sqrt(del);
   area_gr *= sqrt(  (rho2^2 - a^2*del)  );
     
   return area_gr;
}

%% proper area (for a resting observer)
define cal_prop_area(r,dr,a){
   return cal_a_gr(r,a)*dr;
}
    