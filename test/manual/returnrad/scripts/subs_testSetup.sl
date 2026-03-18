require("isisscripts");

variable EXIT_SUCCESS = 1;
variable EXIT_FAILURE = 0;

variable dir_plots = "plots/";
variable dir_refdata = "refdata/";

define is_approx_equal(arr1, arr2, prec){
   variable thres = qualifier("thres",1e-16);
   variable ind1 = where(arr1 > thres);
   variable ind2 = where(arr2 > thres);
   variable n = 1.0*length(arr1);

   variable mean_dev = mean(arr1[ind1]/arr2[ind1]);
   return (abs(mean_dev-1.0)/n < prec )  ? 1 : 0;
}