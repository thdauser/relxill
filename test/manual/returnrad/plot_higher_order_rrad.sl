require("isisscripts");
require("create_output_files.sl");
require("scripts/subs_returnrad.sl");
require("scripts/subs_filenames.sl");

%%%%%
variable Tin = 1.0;
variable _ff = "relxillBB_NS";
%%%%
variable PL_W = 17;
variable PL_H = 10;
variable cols = CB_COLOR_SCHEME2_NB;
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

   variable integ_eflu_prim = sum(val_prim*(ehi-elo)/(0.5*elo+ehi));
   variable integ_eflu_refl = sum(val_refl*(ehi-elo)/(0.5*elo+ehi));;
   
   variable emis = read_ion_gradient(fstruct.disk.rrad);
   
   return struct{elo=elo, ehi=ehi, flux_prim=val_prim, flux_refl=val_refl, emis=emis,
   integ_flu=struct{prim=integ_eflu_prim, refl=integ_eflu_refl}};
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

define print_integ_flu(d){
   variable ii, n = length(d);

   vmessage(" Order #0 energy flux [keV/cm^2/s] : %.2e ", d[0].integ_flu.prim/d[0].integ_flu.prim);
   variable diff;
   _for ii(0, n-1){
      if (ii>0){
	 diff = (d[ii].integ_flu.refl-d[ii-1].integ_flu.refl)/d[ii-1].integ_flu.refl;
      } else {
	 diff = d[ii].integ_flu.refl/d[ii].integ_flu.prim;
      }
      vmessage(" Order #%d energy_flux [keV/cm^2/s] : %.2e (relativ to last order %.1f%%) ",
      ii+1, d[ii].integ_flu.refl/d[0].integ_flu.prim, diff*100);
   }
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
   variable ener_flu_relative = Double_Type[norder];
   _for ii(0, norder-1){
      if (qualifier_exists("order")){
	 lab[ii] = sprintf("order = %d",d[ii].order);
      } else {
	 lab[ii] = sprintf("a = %.3f, incl=%.0f",d[ii].spin, d[ii].incl);
      }

      if (qualifier_exists("order")){
	 variable diff;
	 if (ii==0){
	    diff = d[ii].flux_refl;
	 } else {
	    diff = d[ii].flux_refl-d[ii-1].flux_refl;
	 }
	 
	 pl.plot(emn,diff; color=cols[ii],width=2,line=1);
      }
      pl.plot(emn,d[ii].flux_refl; color=cols[ii],width=2,line=0);
   }

   if (qualifier_exists("order")){
      print_integ_flu(d);
   }

   pl.add_object(xfig_new_legend(lab,cols, 0,2,0.5; labelsize="\small"R),
   0.96,0.96,0.5,0.5;world0);


   return pl;
}

define plot_higher_order_emis(d){
   %%%%%
   variable pl = xfig_plot_new(PL_W, PL_H);

   pl.world(1,50,1e-7,6e-6;loglog);
   
   variable norder = length(d);

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

      pl.plot(rad, emis[ii]; color=cols[ii],width=2);
      pl.plot(rad, diff; color=cols[ii],width=2,line=1);

   }

   pl.xlabel("Radius $[R_\mathrm{g}]$"R);
   pl.ylabel("Emissivity * Radius$^3$");
   
   return pl;
}

%% decide which parameter to plot by checking its length
define get_plot_param(spin, incl){

   variable norder = qualifier("norder",3); 
   variable tprofile = qualifier("tprofile",1);

   variable n;
   if (length(spin)>1){
      n = length(spin);
   } else {
      n = length(incl);
   }
   variable darr = Struct_Type[n];
   variable ii;
   variable d;
   _for ii(0,n-1){
      d =  (
      (length(spin)>1) ?
      get_dat_higher_order(spin[ii], incl, norder; tprofile=tprofile )
      :
      get_dat_higher_order(spin, incl[ii], norder; tprofile=tprofile )
      );

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

variable pl = get_full_plot(spin, incl);
render_with_png(pl, "plots/higher_order_mirrored_zero.pdf");

%% param plots
variable spin_arr = [0.998,0.99,0.95,0.9,0.5,0.0];
variable pl_spin = get_plot_param(spin_arr, incl);
render_with_png(pl_spin, "plots/higher_order_mirrored_spec_spin.pdf");

variable incl_arr = [20.,40.,60.,75.];
render_with_png(
get_plot_param(0.998, incl_arr),
"plots/higher_order_mirrored_spec_incl_a0998.pdf"
);

variable incl_arr = [20.,40.,60.,75.];
render_with_png(
get_plot_param(0.9, incl_arr),
"plots/higher_order_mirrored_spec_incl_a0900.pdf"
);

%% order plots
render_with_png(
get_full_plot(spin, incl;tprofile=tprofile),
"plots/higher_order_mirrored_spec_emis_zero.pdf"
);

variable tprofile = 2;
render_with_png(
get_full_plot(spin, incl; tprofile=tprofile),
"plots/higher_order_mirrored_spec_emis_diskbb.pdf"
);


