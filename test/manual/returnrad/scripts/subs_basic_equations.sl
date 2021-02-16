require("isisscripts");


define cal_a_gr(r,a){
   variable rho2 = r^2 + a^2;  %% this is rho^2
   variable del = r^2 - 2*r + a^2;

   variable area_gr = 2*PI / sqrt(del);
   area_gr *= sqrt(  (rho2^2 + 2*a^2*r)  );
     
   return area_gr;
}


define calc_proper_area_ring(rlo,rhi,a){
   
   variable rmean = 0.5*(rlo+rhi);
   
   variable area_gr = cal_a_gr(rmean,a)*(rhi-rlo);
   return area_gr;
}

define calc_proper_area_ring_relxill(rlo,rhi,a){
   variable factor_relxill = 1./4;
   return calc_proper_area_ring(rlo,rhi,a) * factor_relxill;
   
}