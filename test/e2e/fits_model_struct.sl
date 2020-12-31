require("isisscripts");

define std_eval_model(){
   variable std_elo, std_ehi;
   (std_elo, std_ehi) = log_grid(0.1,1000,2000);
   variable elo = qualifier("elo",std_elo);
   variable ehi = qualifier("ehi",std_ehi);
   
   variable value = eval_fun_keV(elo,ehi);
   return std_elo, std_ehi, value;
}

define fits_write_model_struct(fname){
   
   variable elo, ehi, value;
   (elo, ehi, value) = std_eval_model(;;__qualifiers() );
   
   variable ff = get_fit_fun();
   variable params = reduce_struct(merge_struct_arrays(get_params()),["name","value","freeze","min","max"];extract);
   variable dat = struct{bin_lo=elo,bin_hi=ehi,value=value};

   if (qualifier_exists("verbose")){
         vmessage(" - creating: %s ", fname);
   }
   
   variable fp = fits_open_file(fname,"c");
   fits_write_binary_table(fp,"MODEL",params);   
   fits_update_key(fp,"fitfun",ff);
   
   fits_write_binary_table(fp,"DATA",dat);
   
   fits_close_file(fp);       
}

define fits_read_model_struct(fname){
   
   variable mo = fits_read_table(fname+"[MODEL]");
   variable ff = fits_read_key(fname+"[MODEL]","fitfun");
   variable dat = fits_read_table(fname+"[DATA]");
   
   
   
   variable ii, n = length(mo.name);
   
   fit_fun(ff);
   _for ii(0, n-1, 1){
      set_par(mo.name[ii],mo.value[ii],mo.freeze[ii],mo.min[ii],mo.max[ii]);      
   }
     
   return dat;
}
