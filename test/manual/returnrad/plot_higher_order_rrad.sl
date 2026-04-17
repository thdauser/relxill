require("isisscripts");
require("create_output_files.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

%%%%%
variable Tin = 1.0;
variable _ff = "relxillBB";
%%%%
variable PL_W = 17;
variable PL_H = 10;
%%%

variable fstruct = getFilenameStruct();

define eval_model_all(){
   variable elo, ehi;
   (elo, ehi) = log_grid(0.3,150,500);
   variable fac = (0.5*(elo+ehi))^2/(ehi-elo);

   set_par("*.boost",0);
   variable val_prim = eval_fun_keV(elo,ehi)*fac;

   set_par("*.boost",-1,0,-1,10);
   variable val_refl = eval_fun_keV(elo,ehi)*fac;

   vmessage(" sum prim: %g  \t sum refl: %g ",sum(val_prim*(ehi-elo)),sum(val_refl*(ehi-elo)));

   variable emis = read_ion_gradient(fstruct.disk.rrad);
   
   return struct{elo=elo, ehi=ehi, flux_prim=val_prim, flux_refl=val_refl, emis=emis};
}

define get_spec_order(spin, _incl, norder){

   variable tprofile = qualifier("tprofile",1);
   
   %% does not matter for the mirrored reflection
   variable logxi = 2.0;
   variable logn = 17;

   
   variable params = get_param_struct(spin, Tin, logxi, logn; incl=_incl);
   clean_output_files();
   load_std_model(_ff, params);

   putenv("RELXILL_RRAD_BBODY_TPROFILE=$tprofile"$);
   putenv("RELXILL_RRAD_BBODY_MIRRORED=1");
   putenv("RELXILL_RRAD_BBODY_ORDER=$norder"$);
   putenv("RELXILL_WRITE_FILES=1");
   variable d = eval_model_all();

   return struct_combine(d, struct{order=norder, spin=spin, incl=_incl});
}

define get_dat_higher_order(spin, _incl, norder){

   variable d = Struct_Type[norder];
   variable order = [1:norder];
   variable ii;
   _for ii(0, norder-1){
      d[ii] = get_spec_order(spin, _incl, order[ii];;__qualifiers());
   }

   return d;
}

define plot_higher_order_spec(d){
   %%%%%
   variable pl = xfig_plot_new(PL_W, PL_H);

   variable norder = length(d);
   
   pl.world(d[0].elo[0], d[0].ehi[-1],1e-2,200;loglog);
   
   pl.xlabel("Energy [keV]");
   pl.ylabel("Energy Flux $\nu F_\nu$  [keV$^2$/ cm$^2$\,s\,keV]"R);


   % make life easier
   variable emn = 0.5*(d[0].elo+d[0].ehi);

   pl.plot(emn,d[0].flux_prim; color="black",width=2,line=0);

   variable ii;
   variable lab = String_Type[norder];
   _for ii(0, norder-1){
      if (qualifier_exists("order")){
	 lab[ii] = sprintf("order = %d",d[ii].order);
      } else {
	 lab[ii] = sprintf("a = %.3f, incl=%.0f",d[ii].spin, d[ii].incl);
      }

      if (not qualifier_exists("nodiff")){
	 variable diff;
	 if (ii==0){
	    diff = d[ii].flux_refl;
	 } else {
	    diff = d[ii].flux_refl-d[ii-1].flux_refl;
	 }
	 pl.plot(emn,diff; color=CB_COLOR_SCHEME2_NB[ii],width=2,line=1);
      }
      pl.plot(emn,d[ii].flux_refl; color=CB_COLOR_SCHEME2_NB[ii],width=2,line=0);
   }

   pl.add_object(xfig_new_legend(lab,CB_COLOR_SCHEME_NB, 0,2,0.5; labelsize="\small"R),
   0.96,0.96,0.5,0.5;world0);


   return pl;
}

define plot_higher_order_emis(d){
   %%%%%
   variable pl = xfig_plot_new(PL_W, PL_H);

   pl.world(1,50,1e-7,6e-6;loglog);
   
   variable norder = length(d);

   print(d[0].emis);
   
   variable ii;
   variable emis = Array_Type[norder];
   _for ii(0, norder-1){
      variable rad = 0.5*(d[ii].emis.rlo+d[ii].emis.rhi);
      emis[ii] = d[ii].emis.emis*rad^3;

      variable diff;
      if (ii==0){
	 diff = emis[ii];
      } else {
	 diff = emis[ii] - emis[ii-1];
      }

      pl.plot(rad, emis[ii]; color=CB_COLOR_SCHEME2_NB[ii],width=2);
      pl.plot(rad, diff; color=CB_COLOR_SCHEME2_NB[ii],width=2,line=1);

   }

   pl.xlabel("Radius $[R_\mathrm{g}]$"R);
   pl.ylabel("Emissivity * Radius$^3$");
   
   return pl;
}

define get_plot_spin(spin, incl){

   variable norder = qualifier("norder",3); 
   variable tprofile = qualifier("tprofile",1);

   variable n = length(spin);
   variable darr = Struct_Type[n];
   variable ii;
   _for ii(0,n-1){
      variable d = get_dat_higher_order(spin[ii], incl, norder; tprofile=tprofile );
      darr[ii] = d[-1]; %% only pick the highest order
   }
   
   return plot_higher_order_spec(darr ;nodiff);
}

define get_full_plot(spin, incl){

   variable norder = qualifier("norder",3); 
   variable tprofile = qualifier("tprofile",1);
   
   variable d = get_dat_higher_order(spin, incl, norder; tprofile=tprofile );
   
   variable pl_spec = plot_higher_order_spec(d; order);
   variable pl_emis = plot_higher_order_emis(d);
   
   pl_spec.title(sprintf("spin=%.3f, incl=%.0f, Tin=%.1f, Temp-Profile: %s",
     spin, incl, Tin, (tprofile==1? "zero-torque" : "diskbb" )));

   return  xfig_new_vbox_compound(pl_spec, pl_emis, 0.5);
}


variable spin = 0.998;
variable incl = 40.0;
variable tprofile = 1;
get_full_plot(spin, incl; tprofile=tprofile).render("plots/higher_order_mirrored_spec_emis_zero.pdf");

variable tprofile = 2;
get_full_plot(spin, incl; tprofile=tprofile).render("plots/higher_order_mirrored_spec_emis_diskbb.pdf");

variable spin_arr = [0.998,0.99,0.95,0.9];
get_plot_spin(spin_arr, incl).render("plots/higher_order_mirrored_spec_spin.pdf");