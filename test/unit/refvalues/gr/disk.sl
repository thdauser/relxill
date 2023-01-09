define ad_E(r,a){
   
   return (r*sqrt(r)-2*sqrt(r)+a) / 
     (r^(0.75) *  sqrt( r*sqrt(r) - 3*sqrt(r) + 2*a)  );
   
}
define ad_Lz(r,a){
   
   return sign(a) * (r^2 - 2*a*sqrt(r) + a^2) / 
     (r^(0.75) *  sqrt( r*sqrt(r) - 3*sqrt(r) + 2*a)  );
   
}

define ad_lambda(r,a){
   
   return sign(a) * (r^2 - 2*a*sqrt(r) + a^2) / 
    (r*sqrt(r)-2*sqrt(r)+a) ;   
}



%%% --------------------------- %%%

define qrad_disk(r,rin){
      return (1.-1./sqrt(r/rin)) / r^3;
}


