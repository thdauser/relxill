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

    Copyright 2019 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "test_relxill.h"

void set_std_param_xillver(double *inp_par) {
    inp_par[0] = 2.1;    // Gamma
    inp_par[1] = 1.0;    // Afe
    inp_par[2] = 300.0;  // Ecut
    inp_par[3] = 0.0;    // logxi
    inp_par[4] = 0.;     // redshift
    inp_par[5] = 45.0;   // inclination
    inp_par[6] = -1.0;   // refl. frac.
}

void set_std_param_relline(double *inp_par) {
    inp_par[0] = 6.4;
    inp_par[1] = 3.0;
    inp_par[2] = 3.0;
    inp_par[3] = 15.0;
    inp_par[4] = 0.998;
    inp_par[5] = 30.0;
    inp_par[6] = -1.1;
    inp_par[7] = 1000.;
    inp_par[8] = 0.0;
    inp_par[9] = 1.0;
}

void set_std_param_relconv(double *inp_par) {
    inp_par[0] = 3.0;
    inp_par[1] = 3.0;
    inp_par[2] = 15.0;
    inp_par[3] = 0.998;
    inp_par[4] = 30.0;
    inp_par[5] = -1.;
    inp_par[6] = 400.;
    inp_par[7] = 0.0;
}

void set_std_param_relconvlp(double *inp_par) {
    inp_par[0] = 3.0;
    inp_par[1] = 0.998;
    inp_par[2] = 30.0;
    inp_par[3] = -1.;
    inp_par[4] = 400.;
    inp_par[5] = 0.0;
    inp_par[6] = 2.0;
}

void set_std_param_relxill(double *inp_par) {
    inp_par[0] = 3.0;
    inp_par[1] = 3.0;
    inp_par[2] = 15.0;
    inp_par[3] = 0.998;
    inp_par[4] = 60.0;
    inp_par[5] = -1.0;
    inp_par[6] = 400.;
    inp_par[7] = 0.0;   // redshift
    inp_par[8] = 2.1;   // pl Index
    inp_par[9] = 0.0;   // logxi
    inp_par[10] = 1.0;   // Afe
    inp_par[11] = 300.0; // Ecut
    inp_par[12] = 3.0;   // refl_frac
}

void set_std_param_relxill_nthcomp(double *inp_par) {
    inp_par[0] = 3.0;
    inp_par[1] = 3.0;
    inp_par[2] = 15.0;
    inp_par[3] = 0.998;
    inp_par[4] = 60.0;
    inp_par[5] = -1.0;
    inp_par[6] = 400.;
    inp_par[7] = 0.0;   // redshift
    inp_par[8] = 2.1;   // pl Index
    inp_par[9] = 0.0;   // logxi
    inp_par[10] = 1.0;   // Afe
    inp_par[11] = 100.0; // kTe
    inp_par[12] = 0.0;   // refl_frac
}

void set_std_param_relxilldens(double *inp_par) {
    inp_par[0] = 3.0;
    inp_par[1] = 3.0;
    inp_par[2] = 15.0;
    inp_par[3] = 0.998;
    inp_par[4] = 60.0;
    inp_par[5] = -1.0;
    inp_par[6] = 400.;
    inp_par[7] = 0.0;   // redshift
    inp_par[8] = 2.1;   // pl Index
    inp_par[9] = 0.0;   // logxi
    inp_par[10] = 1.0;   // Afe
    inp_par[11] = 15.0; // logN
    inp_par[12] = 3.0;   // refl_frac
}

void set_std_param_relxilllp(double *inp_par) {
    inp_par[0] = -1.1;   // height
    inp_par[1] = 0.9798; // a
    inp_par[2] = 60.0;  // incl
    inp_par[3] = -1.0;  // rin
    inp_par[4] = 1000.;  // rout
    inp_par[5] = 0.0;    // redshift
    inp_par[6] = 2.1;   // pl Index
    inp_par[7] = 0.0;   // logxi
    inp_par[8] = 1.0;   // Afe
    inp_par[9] = 300.0; // Ecut
    inp_par[10] = 3.0;   // refl_frac
    inp_par[11] = 0.0;   // fixReflFrac
}

void set_std_param_relxilllp_nthcomp(double *inp_par) {
    inp_par[0] = -1.1;   // height
    inp_par[1] = 0.9798; // a
    inp_par[2] = 60.0;  // incl
    inp_par[3] = -1.0;  // rin
    inp_par[4] = 400.;  // rout
    inp_par[5] = 0.0;    // redshift
    inp_par[6] = 2.1;   // pl Index
    inp_par[7] = 0.0;   // logxi
    inp_par[8] = 1.0;   // Afe
    inp_par[9] = 100.0; // kTe
    inp_par[10] = 3.0;   // refl_frac
    inp_par[11] = 0.0;   // fixReflFrac
}

void set_std_param_relxilllpdens(double *inp_par) {
    inp_par[0] = 3.0;   // height
    inp_par[1] = 0.998; // a
    inp_par[2] = 60.0;  // incl
    inp_par[3] = -1.0;  // rin
    inp_par[4] = 400.;  // rout
    inp_par[5] = 0.0;    // redshift
    inp_par[6] = 2.1;   // pl Index
    inp_par[7] = 0.0;   // logxi
    inp_par[8] = 1.0;   // Afe
    inp_par[9] = 16.0; // logN
    inp_par[10] = -1.0;   // refl_frac
    inp_par[11] = 1.0;   // fixReflFrac
}


void set_std_param_relline_lp(double *inp_par) {
    inp_par[0] = 1.0;
    inp_par[1] = 3.0;
    inp_par[2] = 0.998;
    inp_par[3] = 30.0;
    inp_par[4] = -1.;
    inp_par[5] = 400.;
    inp_par[6] = 0.0;  // redshift
    inp_par[7] = 0.0;
    inp_par[8] = 2.0;  // gamma
}

void set_std_param_relxilllpion(double *inp_par) {
    inp_par[0] = 6;   // height
    inp_par[1] = 0.9; // a
    inp_par[2] = 30.0;  // incl
    inp_par[3] = -1.0;  // rin
    inp_par[4] = 1000.;  // rout
    inp_par[5] = 0.0;    // redshift
    inp_par[6] = 2.0;   // pl Index
    inp_par[7] = 3.00;   // logxi
    inp_par[8] = 1.0;   // Afe
    inp_par[9] = 300.0; // Ecut
    inp_par[10] = 0.0;     // beta
    inp_par[11] = 0;     // ion_grad_type
    inp_par[12] = 0.0;   // ion_grad_index
    inp_par[13] = 3.0;   // refl_frac
    inp_par[14] = 0.0;   // fixReflFrac
}

void set_std_param_relxilllpion_nthcomp(double *inp_par) {
    inp_par[0] = 6;   // height
    inp_par[1] = 0.9; // a
    inp_par[2] = 30.0;  // incl
    inp_par[3] = -1.0;  // rin
    inp_par[4] = 1000.;  // rout
    inp_par[5] = 0.0;    // redshift
    inp_par[6] = 2.0;   // pl Index
    inp_par[7] = 3.00;   // logxi
    inp_par[8] = 1.0;   // Afe
    inp_par[9] = 100.0; // kTe
    inp_par[10] = 0.0;     // beta
    inp_par[11] = 0;     // ion_grad_type
    inp_par[12] = 0.0;   // ion_grad_index
    inp_par[13] = 3.0;   // refl_frac
    inp_par[14] = 0.0;   // fixReflFrac
}


xillParam *get_std_param_xillver(int *status) {
    int n_param = NUM_PARAM_XILLVER;
    double *inp_par = (double *) malloc(sizeof(double) * n_param);
    CHECK_MALLOC_RET_STATUS(inp_par, status, NULL)
    set_std_param_xillver(inp_par);
    return init_par_xillver(inp_par, n_param, status);
}


/** standard evaluation of the relline model **/
void std_eval_relline(int *status, int n) {

    printf("\n ==> Evaluating RELLINE MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELLINE;
    double inp_par[NUM_PARAM_RELLINE];
    set_std_param_relline(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 300;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.5, 8.0);

    /* call the relline model */
    double photar[n_ener];
    int ii;
    if (n > 2) {
        for (ii = 0; ii < n; ii++) {
            inp_par[4] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
            tdrelline(ener, n_ener, photar, inp_par, n_param, status);
        }
    } else {
        inp_par[8] = 0.0;
        tdrelline(ener, n_ener, photar, inp_par, n_param, status);
        inp_par[8] = 1.0;
        tdrelline(ener, n_ener, photar, inp_par, n_param, status);
        inp_par[8] = 2.0;
        tdrelline(ener, n_ener, photar, inp_par, n_param, status);
        inp_par[8] = 0.0;
        tdrelline(ener, n_ener, photar, inp_par, n_param, status);
    }
}


/** standard evaluation of the relline model **/
void std_eval_relconv(int *status, int n) {

    printf("\n ==> Evaluating RELCONV MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELCONV;
    double inp_par[NUM_PARAM_RELCONV];
    set_std_param_relconv(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 2000;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.05, 10.0);

    /* call the relline model */
    double photar[n_ener];

    int ii;
    double ener_line = 1.0;
    for (ii = 0; ii < n_ener; ii++) {
        photar[ii] = 0.0;
        if ((ener[ii] < ener_line) && ener[ii + 1] > ener_line) {
            photar[ii] = 1.0; // ener[1]-ener[0];
        }
    }

    // test output
    // save_xillver_spectrum(ener,photar,n_ener,"test_relconv_inp_spectrum.dat");
    tdrelconv(ener, n_ener, photar, inp_par, n_param, status);
    // save_xillver_spectrum(ener,photar,n_ener,"test_relconv_out_spectrum.dat");
}


/** standard evaluation of the relline model **/
void std_eval_relconvlp(int *status, int n) {

    printf("\n ==> Evaluating RELCONV_LP MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELCONVLP;
    double inp_par[NUM_PARAM_RELCONVLP];
    set_std_param_relconvlp(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 2000;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.05, 10.0);

    /* call the relline model */
    double photar[n_ener];

    int ii;
    double ener_line = 1.0;
    for (ii = 0; ii < n_ener; ii++) {
        photar[ii] = 0.0;
        if ((ener[ii] < ener_line) && ener[ii + 1] > ener_line) {
            photar[ii] = 1.0; // ener[1]-ener[0];
        }
    }

    // test output
    // save_xillver_spectrum(ener,photar,n_ener,"test_relconv_inp_spectrum.dat");
    tdrelconvlp(ener, n_ener, photar, inp_par, n_param, status);
    save_xillver_spectrum(ener, photar, n_ener, "test_relconv_out_spectrum.dat");
}


/** standard evaluation of the relxill model **/
void std_eval_relxill(int *status, int n) {

    printf("\n ==> Evaluating RELXILL MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELXILL;
    double inp_par[NUM_PARAM_RELXILL];
    set_std_param_relxill(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 3000;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

    /* call the relline model */
    double photar[n_ener];

    int ii;
    for (ii = 0; ii < n; ii++) {
        if (n > 1) {
            inp_par[3] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
            inp_par[7] = 1.0 * ii / (n - 1) * 0.1 - 0.0;
            inp_par[7] = 0.0;
            inp_par[9] = 1.0 * ii / (n - 1) * 0.1;
            tdrelxill(ener, n_ener, photar, inp_par, n_param, status);
        }
        // printf(" testing a=%.3f , lxi=%.2f \n",inp_par[3],inp_par[9]);
        //printf(" testing z=%.3f \n",inp_par[7]);
    }
    tdrelxill(ener, n_ener, photar, inp_par, n_param, status);

}

/** standard evaluation of the relxill model **/
void std_eval_relxill_nthcomp(int *status, int n) {

    printf("\n ==> Evaluating RELXILL NTHCOMP MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELXILL;
    double inp_par[NUM_PARAM_RELXILL];
    set_std_param_relxill_nthcomp(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 3000;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

    /* call the relline model */
    double photar[n_ener];

    /* int ii;
     for (ii=0; ii<n; ii++){
        if (n>1){
            //			inp_par[3] = 1.0*ii/(n-1)*0.998*2 - 0.998;
            inp_par[7] = 1.0*ii/(n-1)*0.1 - 0.0;
            inp_par[7] = 0.0;
            // inp_par[9] = 1.0*ii/(n-1)*0.1;
        }
        // printf(" testing a=%.3f , lxi=%.2f \n",inp_par[3],inp_par[9]);
        printf(" testing z=%.3f \n",inp_par[7]);
    }*/
    tdrelxill_nthcomp(ener, n_ener, photar, inp_par, n_param, status);
    save_xillver_spectrum(ener, photar, n_ener, "test_relxillCp_spectrum.dat");
}


/** standard evaluation of the relxillD model **/
void std_eval_relxilldens(int *status, int n) {

    printf("\n ==> Evaluating RELXILL_DENS MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELXILL;
    double inp_par[NUM_PARAM_RELXILL];
    set_std_param_relxilldens(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 4000;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

    /* call the relline model */
    double photar[n_ener];
    int ii;
    for (ii = 0; ii < n; ii++) {
        if (n > 1) {
            inp_par[3] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
            inp_par[9] = 1.0 * ii / (n - 1) * 4.7;
        }
        tdrelxilldens(ener, n_ener, photar, inp_par, n_param, status);
    }

}


/** standard evaluation of the relxill model **/
void std_eval_relxilllp(int *status, int n) {

    printf("\n ==> Evaluating RELXILLLP MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELXILLLP;
    double inp_par[NUM_PARAM_RELXILLLP];
    set_std_param_relxilllp(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 100;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

    /* call the relline model */
    double photar[n_ener];
    int ii;


    struct timeval start, end;
    long seconds, useconds;
    gettimeofday(&start, NULL);


    for (ii = 0; ii < n; ii++) {
        if (n > 1) {
            inp_par[1] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
            //		inp_par[7] = 1.0*ii/(n-1)*4.7;
            printf(" relxilllp: testing a=%.3f , lxi=%.2f \n", inp_par[1], inp_par[7]);
        }
        tdrelxilllp(ener, n_ener, photar, inp_par, n_param, status);
    }
    double sum = 0.0;
    int jj;
    for (jj = 0; jj < n_ener; jj++) {
        sum += photar[jj];
    }
    printf(" integ flux = %.2e \n", sum);

    gettimeofday(&end, NULL);


    if (n > 1) {
        seconds = end.tv_sec - start.tv_sec;
        useconds = end.tv_usec - start.tv_usec;

        double mtime = ((seconds) * 1000. + useconds * 0.001) / ((double) n);

        printf("time per relxilllp evaluation: %.1f milli seconds\n", mtime);
    }

    // test output
    // save_xillver_spectrum(ener,photar,n_ener,"test_relxilllp_spectrum.dat");

}


/** standard evaluation of the relxill model **/
void std_eval_relxilllpion(int *status, int n) {

    printf("\n ==> Evaluating RELXILLLPION MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELXILLLPION;
    double inp_par[NUM_PARAM_RELXILLLPION];
    set_std_param_relxilllpion(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 1000;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

    /* call the relline model */
    double photar[n_ener];

    struct timeval start, end;
    long seconds, useconds;
    gettimeofday(&start, NULL);

    int ii;
    for (ii = 0; ii < n; ii++) {
        if (n > 1) {
            inp_par[1] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
            //		inp_par[7] = 1.0*ii/(n-1)*4.7;
            printf(" relxilllpion: testing a=%.3f , lxi=%.2f \n", inp_par[1], inp_par[7]);
        }
        tdrelxilllpion(ener, n_ener, photar, inp_par, n_param, status);
    }

    gettimeofday(&end, NULL);


    if (n > 1) {
        seconds = end.tv_sec - start.tv_sec;
        useconds = end.tv_usec - start.tv_usec;

        double mtime = ((seconds) * 1000 + useconds * 0.001) / ((double) n);

        printf("time per relxilllpion evaluation: %.1f milli seconds\n", mtime);
    }

}


/** standard evaluation of the relxill model **/
void std_eval_relxilllp_nthcomp(int *status, int n) {

    printf("\n ==> Evaluating RELXILLLP NTHCOMP MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELXILLLP;
    double inp_par[NUM_PARAM_RELXILLLP];
    set_std_param_relxilllp_nthcomp(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 100;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

    /* call the relline model */
    double photar[n_ener];
    int ii;

    for (ii = 0; ii < n; ii++) {
        if (n > 1) {
            inp_par[1] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
            //		inp_par[7] = 1.0*ii/(n-1)*4.7;
            printf(" relxilllp (nthcomp): testing a=%.3f , lxi=%.2f \n", inp_par[1], inp_par[7]);
        }
        tdrelxilllp_nthcomp(ener, n_ener, photar, inp_par, n_param, status);

    }

    // save_xillver_spectrum(ener,photar,n_ener,"test_relxilllp_spectrum.dat");

}


/** standard evaluation of the relxill model **/
void std_eval_relxilllpion_nthcomp(int *status, int n) {

    printf("\n ==> Evaluating RELXILLLP ION NTHCOMP MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELXILLLPION;
    double inp_par[NUM_PARAM_RELXILLLPION];
    set_std_param_relxilllpion_nthcomp(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 100;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

    /* call the relline model */
    double photar[n_ener];
    int ii;

    for (ii = 0; ii < n; ii++) {
        if (n > 1) {
            inp_par[1] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
            printf(" relxilllpion (nthcomp): testing a=%.3f , lxi=%.2f \n", inp_par[1], inp_par[7]);
        }
        tdrelxilllpion_nthcomp(ener, n_ener, photar, inp_par, n_param, status);
        printf(" -> %e \n", photar[0]);
    }

}


/** standard evaluation of the relxill model **/
void std_eval_relxilllpdens(int *status, int n) {

    printf("\n ==> Evaluating RELXILLLP_DENS MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELXILLLP;
    double inp_par[NUM_PARAM_RELXILLLP];
    set_std_param_relxilllpdens(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 100;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.1, 1000.0);

    /* call the relline model */
    double photar[n_ener];
    int ii;


    struct timeval start, end;
    long seconds, useconds;
    gettimeofday(&start, NULL);

    for (ii = 0; ii < n; ii++) {
        if (n > 1) {
            inp_par[1] = 1.0 * ii / (n - 1) * 0.998 * 2 - 0.998;
            inp_par[7] = 1.0 * ii / (n - 1) * 4.7;
            printf(" relxilllp: testing a=%.3f , lxi=%.2f \n", inp_par[1], inp_par[7]);
        }
        tdrelxilllpdens(ener, n_ener, photar, inp_par, n_param, status);
    }

    gettimeofday(&end, NULL);

    if (n > 1) {
        seconds = end.tv_sec - start.tv_sec;
        useconds = end.tv_usec - start.tv_usec;

        double mtime = ((seconds) * 1000 + useconds * 0.001) / ((double) n);

        printf("time per relxilllp evaluation: %.1f milli seconds\n", mtime);
    }

}

/** standard evaluation of the relline model **/
void std_eval_relline_lp(int *status, int n) {

    printf("\n ==> Evaluating RELLINE LP MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_RELLINELP;
    double inp_par_lp[n_param];
    set_std_param_relline_lp(inp_par_lp);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 600;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.1, 1.5);

    /* call the relline model */
    double photar[n_ener];
    int ii;
    for (ii = 0; ii < n; ii++) {
//		tdrellinelp(ener,n_ener,photar,inp_par_lp,n_param,status);
        inp_par_lp[2] = 0.998;
        tdrellinelp(ener, n_ener, photar, inp_par_lp, n_param, status);
        tdrellinelp(ener, n_ener, photar, inp_par_lp, n_param, status);
        inp_par_lp[2] = 0.99;
        tdrellinelp(ener, n_ener, photar, inp_par_lp, n_param, status);
        inp_par_lp[2] = 0.998;
        tdrellinelp(ener, n_ener, photar, inp_par_lp, n_param, status);
        CHECK_STATUS_VOID(*status)
    }


}

/** standard evaluation of the relline model **/
void std_eval_xillver(int *status, int n) {

    printf("\n ==> Evaluating XILLVER MODEL \n");
    /* set the parameters */
    int n_param = NUM_PARAM_XILLVER;
    double inp_par[n_param];
    set_std_param_xillver(inp_par);
    CHECK_STATUS_VOID(*status)

    /* create an energy grid */
    int n_ener = 2000;
    double ener[n_ener + 1];
    get_log_grid(ener, n_ener + 1, 0.08, 500.0);


    /* call the relline model */
    double photar[n_ener];
    int ii;
    for (ii = 0; ii < n; ii++) {
        tdxillver(ener, n_ener, photar, inp_par, n_param, status);
        CHECK_STATUS_VOID(*status)
    }

    // test output
    // save_xillver_spectrum(ener,photar,n_ener,"test_xillver_spectrum.dat");


}

