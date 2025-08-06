variable table_name = "XrB_mu1Mdot31Rj25Risco3Omega2Ms11b1p1_JED.fits";

variable global_fitfun = "jedsad+relxill_jedsad";

variable dir_tables = getenv("RELXILL_TABLE_PATH");

define load_models(){
   
   load_xspec_local_models("./build/");
   __set_hard_limits ("relxill_jedsad","Rin", -100,100; parm=[2,2,6]);
   
   add_atable_model(dir_tables+"/jedsad/$table_name"$, "jedsad");
   
   fit_fun(global_fitfun);
}


define eval_ff(){
   variable elo, ehi;
   (elo, ehi) = log_grid(0.2,200,300);

   variable val = eval_fun_keV(elo,ehi)/(ehi-elo)*(0.5*(elo+ehi))^2;

   if (qualifier_exists("norm")){
      variable norm_energy = qualifier("norm_energy",1.0);
      variable norm_value = qualifier("norm",1.0);
      variable ind = wherefirst(elo>norm_energy);
      val *= norm_value / val[ind];      
   }
   
   return struct{lo=elo, hi=ehi, emean=0.5*(elo+ehi), val=val};
}

define setup_model(){

   fit_fun(global_fitfun);
   set_par_fun("relxill_jedsad(1).Rin","jedsad(1).Risco");  % jedsad only starts at 2Rg

   variable tie_params = ["Rj", "Mdot", "Omega", "Ms"];
   variable p;
   variable param;
   foreach p(tie_params){
      variable pname_rel = "relxill_jedsad(1).$p"$;
      variable pname_jed = "jedsad(1).$p"$;
      param = get_params(pname_rel)[0];
      set_par(pname_jed,param.value, param.freeze, param.min, param.max);
      set_par_fun(pname_rel, pname_jed);
   }
   
%   set_par_fun("jedsad(1).Risco", "abs(relxill_jedsad(1).Rin)*kerr_rms(relxill_jedsad(1).a)");   
}
