require("isisscripts");

variable EXIT_SUCCESS = 1;
variable EXIT_FAILURE = 0;

variable dir_plots = "plots/";
variable dir_refdata = "refdata/";

define is_approx_equal(arr, arr_ref, prec){
   variable thres = qualifier("thres",1e-3);
   variable ind = where(arr_ref/max(arr_ref) > thres);
   
   variable mean_dev = mean(arr[ind]/arr_ref[ind]);
   return (abs(mean_dev-1.0) < prec )  ? 1 : 0;
}