% -*- mode: slang; mode: fold -*-
require("gr");
require("relxill");

%% require("physical_disk_irradiation/integrated_emissivity_flux.sl");

variable CONST_cm2pc = 3.2407792700054E-19;
variable CONST_rgIncm = 1.4822e5;  %% G / c^2 * (M/Msolar)

define convertFobs2Lobs(flux0, D_kpc){ %{{{
   
   %% unit flux0 = erg/cm²/s
   
   variable D_cm = D_kpc * 1e3/ CONST_cm2pc;   
   
   variable lobs = 4 * PI * D_cm^2 * flux0;
   
   return lobs;
}
%}}}

define calcFluxRin(Linc, h_rg, Rin_rg, M){ %{{{
   
   variable dRay_rg = sqrt(h_rg^2 + Rin_rg^2);
   variable dRay_cm = dRay_rg * CONST_rgIncm * M;
   
   variable flux = Linc / (4 * PI * dRay_cm^2);
   
   return flux;
}
%}}}

define calcFluxBoostPrimary(height, spin, gamma){ %{{{
   
   variable ii, n = length(height);
   variable fluxBoost = Double_Type[n];
   
   _for ii(0, n-1){
      fluxBoost[ii] = gr->lp_fluxBoost_obs2primary(height[ii], spin, gamma[ii]);
   }
   
   return fluxBoost;
}
%}}}

define calcLsource(height, spin, fboostPrim, rin, rout, lobs){ %{{{
   
   variable ii, n = length(height);
   variable reflFrac = reflection_fraction_relxill(spin, height, rin, rout;struct);
      
   return  lobs * fboostPrim / (reflFrac.f_inf / 0.5);
}
%}}}





define getEmisPhysnorm(inpDat, setup){ %{{{

   variable ii, n = length(inpDat.height);
   
   variable emisProf = Struct_Type[n];
   
   variable emisNonrelat = Array_Type[n];
   
   _for ii(0, n-1){ 
      emisProf[ii] = getEmisProfile(inpDat.height[ii],setup.a,inpDat.gamma[ii]; Rin=setup.Rin);

      emisNonrelat[ii] = gr->lp_emisProfile_nonrelat(inpDat.height[ii],emisProf[ii].rad);

      variable normFacFlatRin = inpDat.fluxFlatRin[ii] / emisNonrelat[ii][0]; 
	  
      emisProf[ii].emis *= normFacFlatRin*inpDat.reflFrac[ii]*
	inpDat.fluxBoostPrimary[ii]; %% convert from Lobs to Lsrc

      emisNonrelat[ii] *= normFacFlatRin;
      
   }
   

   return emisProf, emisNonrelat;
}
%}}}


%!%+
%\function{calc_physical_irradiation(setup, Lobs) }
%
%!%-
define calc_physical_irradiation(setup, Lobs, Msol){ %{{{

   variable spin = setup.a;
   variable height = setup.h;
   variable rin = setup.Rin;
   variable rout = setup.Rout;
   variable gamma = setup.gamma;

   variable r_lximax = rin*(11./9.)^2;
   
   
   variable reflFrac =  reflection_fraction_relxill(spin,height,rin,rout; struct);

   variable fluxBoostPrimary =  gr->lp_fluxBoost_obs2primary(height, spin, gamma); 
   
   variable Lsrc =  calcLsource(height, spin, fluxBoostPrimary, rin,rout,Lobs);

   variable fluxFlatRin = calcFluxRin(Lsrc, height, rin, Msol);
   variable fluxFlatRlximax = calcFluxRin(Lsrc, height, r_lximax, Msol);
            
   variable physical_irrad_struct = struct{
      model_params=setup,
      r_lximax=r_lximax,
      Lsrc=Lsrc, 
      fluxBoostPrimary=fluxBoostPrimary,
      fluxFlatRin=fluxFlatRin,
      fluxFlatRlximax=fluxFlatRlximax,
      reflFrac=reflFrac.fR, f_inf=reflFrac.f_inf};
   
   return physical_irrad_struct;
}
%}}}


define getEmisPhysnorm_NEW(model_params, physical_irrad){ %{{{

   
   variable emisProf = getEmisProfile(model_params.h,model_params.a,model_params.gamma; Rin=model_params.Rin);
   
   variable emisNonrelat = gr->lp_emisProfile_nonrelat(model_params.h,emisProf.rad) ;
   variable emisNonrelat_r_lximax = gr->lp_emisProfile_nonrelat(model_params.h,physical_irrad.r_lximax) ;

   variable normFacFlatRin = physical_irrad.fluxFlatRin / emisNonrelat[0];
   variable normFacFlatRlximax = physical_irrad.fluxFlatRlximax / emisNonrelat_r_lximax;
  
%   vmessage(" flux_r_lximax=%.4e, emis_newton_r_lximax=%.4e, norm_fac=%.4e",
%	    physical_irrad.fluxFlatRlximax,emisNonrelat_r_lximax,normFacFlatRlximax);

   %% we now use R_lxi_max here (like in the relxilllpAlpha model)
%%   emisProf.emis *= normFacFlatRin*physical_irrad.reflFrac;   
%%   emisNonrelat *= normFacFlatRin; %%/physical_irrad.fluxBoostPrimary; %% to take into account that we already included the boost in the calculation of Lsource

   emisProf.emis *= normFacFlatRlximax*physical_irrad.reflFrac;   
   emisNonrelat *= normFacFlatRlximax; %%/physical_irrad.fluxBoostPrimary; %% to take into account that we already included the boost in the calculation of Lsource

 
   return emisProf, emisNonrelat;
}
%}}}




%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%%%%


