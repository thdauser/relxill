/*
   This file is part of the RELXILL model code.

   RELXILL is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   RELXILL is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.
   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

    Copyright 2016 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "relutility.h"



// TODO: CHECK IF WE NEED "INLINE" HERE for maximal speed
/** linear interpolation in 1 dimension **/
double interp_lin_1d(double ifac_r, double rlo, double rhi){
	return ifac_r*rhi + (1.0-ifac_r)*rlo;
}

double interp_log_1d(double ifac_r, double rlo, double rhi){
	return exp(ifac_r*log(rhi) + (1.0-ifac_r)*log(rlo));
}

// TODO: CHECK IF WE NEED "INLINE" HERE for maximal speed
/** linear interpolation in 2 dimensions **/
double interp_lin_2d(double ifac1, double ifac2, double r11, double r12, double r21, double r22){
	return (1.0 - ifac1) * (1.0 - ifac2) * r11 +
		   (ifac1)       * (1.0 - ifac2) * r12 +
		   (1.0 - ifac1) * (ifac2)       * r21 +
		   (ifac1)       * (ifac2)       * r22;
}

double interp_lin_2d_float(double ifac1, double ifac2, float r11, float r12, float r21, float r22){
	return (1.0 - ifac1) * (1.0 - ifac2) * r11 +
		   (ifac1)       * (1.0 - ifac2) * r12 +
		   (1.0 - ifac1) * (ifac2)       * r21 +
		   (ifac1)       * (ifac2)       * r22;
}

void relxill_error(const char* const func, const char* const msg, int* status){
	*status = EXIT_FAILURE;
	printf(" *** error in relxill (%s): %s!\n", func, msg);
}

void check_relxill_error(const char* const func, const char* const msg, int* status){
	if (*status!=EXIT_SUCCESS){
		*status = EXIT_FAILURE;
		printf(" *** error in relxill (%s): %s!\n", func, msg);
	}
}

void get_version_number(char** vstr, int* status){
	if (asprintf(vstr, "%i.%i.%i", version_major, version_minor, version_build) == -1){
		RELXILL_ERROR("failed to get version number",status);
	}
}

/**  FLOAT search for value "val" in array "arr" (sorted ASCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int binary_search_float(float* arr,int n,float val){

	int klo=0;
	int khi=n-1;
	int k=-1;
	while ( (khi-klo) > 1 ){
		k=(khi+klo)/2;
		if(arr[k]>val){
			khi=k;
		} else {
			klo=k;
		}
	}
	return klo;
}

/**  search for value "val" in array "arr" (sorted ASCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int binary_search(double* arr,int n,double val){

	int klo=0;
	int khi=n-1;
	int k=-1;
	while ( (khi-klo) > 1 ){
		k=(khi+klo)/2;
		if(arr[k]>val){
			khi=k;
		} else {
			klo=k;
		}
	}
	return klo;
}


/**  FLOAT search for value "val" in array "arr" (sorted DESCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int inv_binary_search_float(float* arr,int n,float val){

	int klo=0;
	int khi=n-1;
	int k=-1;
	while ( (khi-klo) > 1 ){
		k=(khi+klo)/2;
		if(arr[k]<val){
			khi=k;
		} else {
			klo=k;
		}
	}
	return klo;
}


/**  search for value "val" in array "arr" (sorted DESCENDING!) with length n and
 	 return bin k for which arr[k]<=val<arr[k+1] **/
int inv_binary_search(double* arr,int n,double val){

	int klo=0;
	int khi=n-1;
	int k=-1;
	while ( (khi-klo) > 1 ){
		k=(khi+klo)/2;
		if(arr[k]<val){
			khi=k;
		} else {
			klo=k;
		}
	}
	return klo;
}

/** test if it is a relxill flavour model **/
int is_relxill_model(int model_type){
	if (model_type < 0){
		return 1;
	} else {
		return 0;
	}
}


/** trapez integration around a single bin **/
double trapez_integ_single(double* re, int ii, int nr){
	if (ii==0){
		return 1.0/re[ii+1] - 1.0/re[ii];
	} else if (ii==nr-1){
        return 1.0/re[ii] - 1.0/re[ii-1];
	} else {
        return 1.0/re[ii+1] - 1.0/re[ii-1];
	}
}

/** convert gstar to energy */
double gstar2ener(double g, double gmin, double gmax, double ener){
	return (g*(gmax-gmin) + gmin)*ener;
}

/** get a radial grid on the accretion disk in order to calculate a relline for each zone **/
double* get_rzone_grid(double rmin, double rmax, int nzones, int* status){

	// rgrid has length nzones+1
	double* rgrid = (double*) malloc( (nzones+1)*sizeof(double));
	CHECK_MALLOC_RET_STATUS(rgrid,status,NULL);

	if (nzones==1){
		rgrid[0] = rmin;
		rgrid[1] = rmax;
	} else {
		get_log_grid(rgrid,nzones+1,rmin,rmax);
//		printf(" *** Warning: A radial zones grid is not yet implemented \n");
	}
	return rgrid;
}

/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_log_grid(double* ener, int n_ener, double emin, double emax){
	for (int ii=0; ii<n_ener; ii++){
		ener[ii] = 1.0*ii / (n_ener-1) * ( log(emax) - log(emin)) + log(emin);
		ener[ii] = exp(ener[ii]);
	}
}


/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_lin_grid(double* ener, int n_ener, double emin, double emax){
	for (int ii=0; ii<n_ener; ii++){
		ener[ii] = 1.0*ii / (n_ener-1) * ( emax - emin) + emin;
	}
}


/* get RMS (ISCO) for the Kerr Case */
double kerr_rms(double a){
	//	 accounts for negative spin
  double sign = 1.0;
  if (a<0) {
     sign = -1.0;
  }

  double Z1 = 1.0+pow(1.0-a*a,1.0/3.0)*(pow(1.0+a,1.0/3.0)+pow(1.0-a,1.0/3.0));
  double Z2=sqrt((3.0*a*a)+(Z1*Z1));

  return 3.0+Z2-sign*sqrt((3.0-Z1)*(3.0+Z1+(2*Z2)));
}


/** calculate the doppler factor for a moving primary source **/
double doppler_factor(double del, double bet) {
	return sqrt(1.0 - bet*bet) / (1.0 + bet*cos(del));
}


/** calculates g = E/E_i in the lamp post geometry (see, e.g., 27 in Dauser et al., 2013, MNRAS) **/
double gi_potential_lp(double r, double a, double h, double bet, double del){

	/** ! calculates g = E/E_i in the lamp post geometry
	  ! (see, e.g., page 48, Diploma Thesis, Thomas Dauser) **/
	double ut_d = ((r*sqrt(r)+a)/(sqrt(r)*sqrt(r*r -3*r + 2*a*sqrt(r))));
	double ut_h = sqrt((h*h + a*a)/(h*h - h + a*a));

	double gi = ut_d/ut_h;

	// check if we need to calculate the additional factor for the velocity
	if (fabs(bet) < 1e-6){
		return gi;
	}

	double gam = 1.0/sqrt(1.0-bet*bet);

	// get the sign for the equation
	double sign = 1.0;
	if (del > M_PI/2) {
		sign = -1.0;
	}

	double delta_eq = h*h - 2*h + a*a;
	double q2 = (pow(sin(del),2))*(  pow((h*h + a*a),2) / delta_eq ) - a*a;

	double beta_fac = sqrt(  pow((h*h + a*a ),2) - delta_eq*(q2 + a*a) );
	beta_fac = gam*(1.0 + sign*beta_fac / (h*h + a*2) *bet);

	return gi / beta_fac;
}



/** print the xillver spectrum   **/
void save_xillver_spectrum(double* ener, double* flu, int n_ener){

	FILE* fp =  fopen ( "test_xillver_spectrum.dat","w+" );
	int ii;
	for (ii=0; ii<n_ener; ii++){
		fprintf(fp, " %e \t %e \t %e \n",ener[ii],ener[ii+1],flu[ii]);
	}
	if (fclose(fp)) exit(1);
}


/* A simple implementation of the FFT taken from http://paulbourke.net/miscellaneous/dft/
   (uses the Radix-2 Cooley-Tukey algorithm)

   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform
*/
void FFT_R2CT(short int dir,long m,double *x,double *y){

   long n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++){
      n *= 2;
   }

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0;
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0;
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1;
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1)
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }

   return;
}


/** rebin spectrum to a given energy grid
 *  length of ener is nbins+1       **/

void rebin_spectrum(double* ener, double* flu, int nbins, double* ener0, double* flu0, int nbins0){

	int ii; int jj;
	int imin = 0;
	int imax = 0;

	for (ii=0; ii<nbins; ii++){

		flu[ii] = 0.0;

		/* check of the bin is outside the given energy range */
		if ( (ener0[0] < ener[ii+1]) && (ener0[nbins0] > ener[ii]) ){

			/* need to make sure we are in the correct bin */
			while ( ener0[imin]<=ener[ii] && imin<=nbins0){
				imin++;
			}
			// need to set it back, as we just crossed to the next bin
			if (imin>0){
				imin--;
			}
			while ( (ener0[imax]<=ener[ii+1] && imax<nbins0)){
				imax++;
			}
			if (imax>0){
				imax--;
			}

			double elo = ener[ii];
			double ehi = ener[ii+1];
			if (elo < ener0[imin]) elo=ener0[imin];
			if (ehi > ener0[imax+1]) ehi=ener0[imax+1];

			if (imax==imin){
				flu[ii] = (ehi-elo) / (ener0[imin+1] - ener0[imin]) * flu0[imin];
			} else {

				double dmin=(ener0[imin+1]-elo)/(ener0[imin+1]-ener0[imin]);
				double dmax=(ehi-ener0[imax])/(ener0[imax+1]-ener0[imax]);

				flu[ii] += flu0[imin]*dmin + flu0[imax]*dmax;

				for (jj=imin+1; jj <= imax-1; jj++) {
					flu[ii] += flu0[jj];
				}

			}

	/**		printf("[%i]  %i-%i  -> ener=[%.3f-%.3f] , choosing bins %.3f and %.3f  => flux=%.2f\n",ii,imin,imax,
								ener[ii],ener[ii+1],ener0[imin],ener0[imax+1],flu[ii]);  **/

		}

	}
}
