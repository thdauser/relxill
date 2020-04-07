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

    Copyright 2020 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "reltable.h"
#include "time.h"

static relDat* new_relDat(int nr, int ng, int* status){
	relDat* dat = (relDat*) malloc (sizeof(relDat));
	CHECK_MALLOC_RET_STATUS(dat,status,dat);

	dat->r = (float*) malloc( sizeof(float) * nr);
	CHECK_MALLOC_RET_STATUS(dat->r,status,dat);
	dat->gmin = (float*) malloc( sizeof(float) * nr);
	CHECK_MALLOC_RET_STATUS(dat->gmin,status,dat);
	dat->gmax = (float*) malloc( sizeof(float) * nr);
	CHECK_MALLOC_RET_STATUS(dat->gmax,status,dat);

	dat->cosne1 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->cosne1,status,dat);
	dat->cosne2 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->cosne2,status,dat);
	dat->trff1 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->trff1,status,dat);
	dat->trff2 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->trff2,status,dat);

	int ii;
	for (ii=0; ii<nr; ii++){
		dat->cosne1[ii] = (float*) malloc( sizeof(float) * ng);
		CHECK_MALLOC_RET_STATUS(dat->cosne1,status,dat);
		dat->cosne2[ii] = (float*) malloc( sizeof(float) * ng);
		CHECK_MALLOC_RET_STATUS(dat->cosne2,status,dat);
		dat->trff1[ii] = (float*) malloc( sizeof(float) * ng);
		CHECK_MALLOC_RET_STATUS(dat->trff1,status,dat);
		dat->trff2[ii] = (float*) malloc( sizeof(float) * ng);
		CHECK_MALLOC_RET_STATUS(dat->trff2,status,dat);
	}
	return dat;
}

/* create a new LP table structure*/
static lpDat* new_lpDat(int n_h, int n_rad, int* status){
	lpDat* dat = (lpDat*) malloc(sizeof(lpDat)*n_h);
	CHECK_MALLOC_RET_STATUS(dat,status,NULL);

	dat->h = (float*) malloc (sizeof(float)*n_h);
	CHECK_MALLOC_RET_STATUS(dat->h,status,dat);
	dat->rad = (float*) malloc (sizeof(float)*n_rad);
	CHECK_MALLOC_RET_STATUS(dat->rad,status,dat);

	dat->intens = (float**) malloc (sizeof(float*)*n_h);
	CHECK_MALLOC_RET_STATUS(dat->intens,status,dat);
	dat->del = (float**) malloc (sizeof(float*)*n_h);
	CHECK_MALLOC_RET_STATUS(dat->del,status,dat);
	dat->del_inc = (float**) malloc (sizeof(float*)*n_h);
	CHECK_MALLOC_RET_STATUS(dat->del_inc,status,dat);

	int ii;
	for (ii=0;ii<n_h;ii++){
		dat->intens[ii] = (float*) malloc (sizeof(float)*n_rad);
		CHECK_MALLOC_RET_STATUS(dat->intens[ii],status,dat);
		dat->del[ii] = (float*) malloc (sizeof(float)*n_rad);
		CHECK_MALLOC_RET_STATUS(dat->del[ii],status,dat);
		dat->del_inc[ii] = (float*) malloc (sizeof(float)*n_rad);
		CHECK_MALLOC_RET_STATUS(dat->del_inc[ii],status,dat);
	}
	return dat;
}

/** get a new and empty rel table (structure will be allocated)  */
relTable* new_relTable(int n_a, int n_mu0, int n_r, int n_g, int* status){
	relTable* tab = (relTable*) malloc (sizeof(relTable));
	CHECK_MALLOC_RET_STATUS(tab,status,tab);

	// we know which dimensions the table should have
	tab->n_a   =  n_a;
	tab->n_mu0 =  n_mu0;
	tab->n_r   =  n_r;
	tab->n_g   =  n_g;

	tab->a = NULL;
	tab->mu0 = NULL;

	tab->arr=NULL;

	tab->arr = (relDat***) malloc (sizeof(relDat**)*tab->n_a);
	CHECK_MALLOC_RET_STATUS(tab->arr,status,tab);

	int ii; int jj;
	for (ii=0; ii<tab->n_a; ii++){

		tab->arr[ii] = (relDat**) malloc (sizeof(relDat*)*tab->n_mu0);
		CHECK_MALLOC_RET_STATUS(tab->arr[ii],status,tab);

		for (jj=0; jj<tab->n_mu0; jj++){
			tab->arr[ii][jj] = NULL;
			CHECK_STATUS_RET(*status,tab);
		}
	}

	return tab;
}

/** read one axis of the rel table from the FITS file   */
static void get_reltable_axis(int nrows, float** val, char* extname, char* colname, fitsfile* fptr, int* status){

	int extver = 0;
	fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %s\n",extname);
		return;
	}

	// get the column id-number
	int colnum;
	if(fits_get_colnum(fptr, CASEINSEN, colname, &colnum, status)) return;

	// get the number of rows
	long n;
	if (fits_get_num_rows(fptr, &n, status)) return;

	if (nrows != n){
		RELXILL_ERROR("wrong dimension of at least one axis given in the rel_table",status);
	}

	// allocate memory for the array
	*val=(float*)malloc(n*sizeof(float));
	CHECK_MALLOC_VOID_STATUS(*val,status);

    int anynul=0;
    double nullval=0.0;
    LONGLONG nelem = (LONGLONG) n;
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nelem ,&nullval,*val, &anynul, status);

	return;
}

/**     */
static void load_single_relDat_2dcol(fitsfile* fptr, float** val,int n1, int n2, int colnum, int* status){

    int anynul=0;
    double nullval=0.0;

    assert(val!=NULL);

    LONGLONG nelem = (LONGLONG) n2;

    int ii;
    for (ii=0; ii<n1;ii++){
        if(fits_read_col(fptr, TFLOAT, colnum, ii+1, 1, nelem ,&nullval,val[ii], &anynul, status)) return;

    }
}

/** load one single data extension from the relline table   */
static relDat* load_single_relDat(fitsfile* fptr, char* extname, int nhdu, int* status){

	// int extver = 0;
	// fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
	int exttype;
	fits_movabs_hdu(fptr,nhdu,&exttype, status);
	if (*status!=EXIT_SUCCESS){
		printf(" *** error moving to extension %s\n",extname);
		return NULL;
	}

	// get the column id-number
	int colnum_r,colnum_gmin,colnum_gmax;
	int colnum_trff1, colnum_trff2;
	int colnum_cosne1, colnum_cosne2;
	if(fits_get_colnum(fptr, CASEINSEN, "r", &colnum_r, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "gmin", &colnum_gmin, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "gmax", &colnum_gmax, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "trff1", &colnum_trff1, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "trff2", &colnum_trff2, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "cosne1", &colnum_cosne1, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "cosne2", &colnum_cosne2, status)) return NULL;

	// check the number of rows (need to coincide with RELTABLE_NR
	long n;
	if (fits_get_num_rows(fptr, &n, status)) return NULL;
	if (n != RELTABLE_NR){
		RELXILL_ERROR("inconsistent number of rows in rel table",status);
		printf("    -> expecting %i, but found %ld in extensions %s",RELTABLE_NR,n,extname);
		return NULL;
	}

	// allocate the memory for the table
	relDat* dat = new_relDat(RELTABLE_NR,RELTABLE_NG,status);
	CHECK_STATUS_RET(*status,NULL);

	// now load the table column by column

    // (1) start with the 1D columns
    int anynul=0;
    double nullval=0.0;
    LONGLONG nelem = (LONGLONG) RELTABLE_NR;
    fits_read_col(fptr, TFLOAT, colnum_r, 1, 1, nelem ,&nullval,dat->r, &anynul, status);
    CHECK_STATUS_RET(*status,dat);
    fits_read_col(fptr, TFLOAT, colnum_gmin, 1, 1, nelem ,&nullval,dat->gmin, &anynul, status);
    CHECK_STATUS_RET(*status,dat);
    fits_read_col(fptr, TFLOAT, colnum_gmax, 1, 1, nelem ,&nullval,dat->gmax, &anynul, status);
    CHECK_STATUS_RET(*status,dat);

    // (2) and finally the 2D columns
    load_single_relDat_2dcol(fptr,dat->trff2,RELTABLE_NR,RELTABLE_NG,colnum_trff2,status);
    CHECK_STATUS_RET(*status,dat);
    load_single_relDat_2dcol(fptr,dat->trff1,RELTABLE_NR,RELTABLE_NG,colnum_trff1,status);
    CHECK_STATUS_RET(*status,dat);
    load_single_relDat_2dcol(fptr,dat->cosne1,RELTABLE_NR,RELTABLE_NG,colnum_cosne1,status);
    CHECK_STATUS_RET(*status,dat);
    load_single_relDat_2dcol(fptr,dat->cosne2,RELTABLE_NR,RELTABLE_NG,colnum_cosne2,status);
    CHECK_STATUS_RET(*status,dat);

    return dat;

}

/** load the complete relline table */
void read_relline_table(char* filename, relTable** inp_tab, int* status){

	relTable* tab = (*inp_tab);
	fitsfile *fptr=NULL;

	char* fullfilename=NULL;
	char* extname=NULL;

	do{ // Errot handling loop
		if (tab != NULL){
			RELXILL_ERROR("relline table already loaded",status);
			break;
		}

		tab = new_relTable(RELTABLE_NA,RELTABLE_NMU0,RELTABLE_NR,RELTABLE_NG,status);
		CHECK_STATUS_BREAK(*status);

		// should be set by previous routine
		assert(tab!=NULL);
		assert(tab->arr!=NULL);

		// get the full filename
		if (asprintf(&fullfilename, "%s/%s", get_relxill_table_path() ,filename) == -1){
			RELXILL_ERROR("failed to construct full path the rel table",status);
			break;
		}

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
			CHECK_RELXILL_ERROR("opening of the rel table failed",status);
			printf("    either the full path given (%s) is wrong \n",fullfilename);
			printf("    or you need to download the table ** %s **  from \n",filename);
			printf("    http://www.sternwarte.uni-erlangen.de/research/relxill/ \n");
			break;
		}

		// first read the axes of the table
		get_reltable_axis(tab->n_a, &(tab->a), "a", "a", fptr, status);
		CHECK_RELXILL_ERROR("reading of spin axis failed",status);

		get_reltable_axis(tab->n_mu0, &(tab->mu0), "mu0", "mu0", fptr, status);
		CHECK_RELXILL_ERROR("reading of mu0 axis failed",status);

		//now load the full table (need to go through all extensions)
		int ii; int jj;
		for (ii=0; ii<tab->n_a; ii++){
			for (jj=0; jj<tab->n_mu0; jj++){

				if (asprintf(&extname, "%i_%i", ii+1,jj+1) == -1){
					RELXILL_ERROR("failed to construct full path the rel table",status);
					break;
				}

				assert(tab->arr[ii][jj]==NULL);
				int nhdu = (ii)*tab->n_mu0+jj+4;
				tab->arr[ii][jj] = load_single_relDat(fptr, extname, nhdu, status);
				free(extname);
				if (*status!=EXIT_SUCCESS){
					RELXILL_ERROR("failed to load data from the rel table into memory",status);
					break;
				}
			}
		}

	} while(0);

	if (*status==EXIT_SUCCESS){
		// assigne the value
		(*inp_tab) = tab;
	} else {
		free_relTable(tab);
	}
	free(fullfilename);

	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}

lpDat* load_single_lpDat(fitsfile* fptr, int n_h, int n_rad, int rownum, int* status){
	lpDat* dat = new_lpDat(n_h,n_rad,status);
	CHECK_MALLOC_RET_STATUS(dat,status,NULL);

	int colnum_r;
	int colnum_hgrid;
	int colnum_h;
	int colnum_del;
	int colnum_del_inc;


	char* colname_h=NULL;
	char* colname_del=NULL;
	char* colname_del_inc=NULL;

	LONGLONG nelem;
    int anynul=0;
    double nullval=0.0;


	if(fits_get_colnum(fptr, CASEINSEN, "r", &colnum_r, status)) return NULL;
	if(fits_get_colnum(fptr, CASEINSEN, "hgrid", &colnum_hgrid, status)) return NULL;


	nelem = (LONGLONG) n_h;
    fits_read_col(fptr, TFLOAT, colnum_hgrid, rownum, 1, nelem ,&nullval,dat->h, &anynul, status);

    nelem = (LONGLONG) n_rad;
    fits_read_col(fptr, TFLOAT, colnum_r, rownum, 1, nelem ,&nullval,dat->rad, &anynul, status);
    CHECK_STATUS_RET(*status,NULL);


	int ii;
    nelem = (LONGLONG) n_rad;
	for (ii=0; ii<n_h; ii++){

		if (asprintf(&colname_h, "h%i", ii+1) == -1){
			RELXILL_ERROR("failed to construct colname of the lp table",status);
			return NULL;
		}
		if(fits_get_colnum(fptr, CASEINSEN, colname_h, &colnum_h, status)) return NULL;
		free(colname_h);

		if (asprintf(&colname_del, "del%i", ii+1) == -1){
			RELXILL_ERROR("failed to construct colname of the lp table",status);
			return NULL;
		}
		if(fits_get_colnum(fptr, CASEINSEN, colname_del, &colnum_del, status)) return NULL;
		free(colname_del);

		if (asprintf(&colname_del_inc, "del_inc%i", ii+1) == -1){
			RELXILL_ERROR("failed to construct colname of the lp table",status);
			return NULL;
		}
		if(fits_get_colnum(fptr, CASEINSEN, colname_del_inc, &colnum_del_inc, status)) return NULL;
		free(colname_del_inc);

	    fits_read_col(fptr, TFLOAT, colnum_h, rownum, 1, nelem ,&nullval,dat->intens[ii], &anynul, status);
	    fits_read_col(fptr, TFLOAT, colnum_del, rownum, 1, nelem ,&nullval,dat->del[ii], &anynul, status);
	    fits_read_col(fptr, TFLOAT, colnum_del_inc, rownum, 1, nelem ,&nullval,dat->del_inc[ii], &anynul, status);

	    CHECK_STATUS_RET(*status,NULL);
	}


	return dat;
}


/** load the complete relline table */
void read_lp_table(char* filename, lpTable** inp_tab, int* status){

	lpTable* tab = (*inp_tab);
	fitsfile *fptr=NULL;

	char* fullfilename=NULL;

	do{ // Errot handling loop
		if (tab != NULL){
			RELXILL_ERROR("relline LP table already loaded",status);
			break;
		}

		tab = new_lpTable(LPTABLE_NA,LPTABLE_NH,LPTABLE_NR,status);
		CHECK_STATUS_BREAK(*status);

		// should be set by previous routine
		assert(tab!=NULL);
		assert(tab->dat!=NULL);

		// get the full filename
		if (asprintf(&fullfilename, "%s/%s", get_relxill_table_path(),filename) == -1){
			RELXILL_ERROR("failed to construct full path the lp table",status);
			break;
		}

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
			CHECK_RELXILL_ERROR("opening of the lp table failed",status);
			printf("    full path given: %s \n",fullfilename);
			break;
		}

		// first read the spin axes of the table AND move to the correct extension
		get_reltable_axis(tab->n_a, &(tab->a), "I_h", "a", fptr, status);
		CHECK_RELXILL_ERROR("reading of spin axis failed",status);


		//now load the full table (need to go through all extensions)
		int rownum=-1;
		int ii;
		for (ii=0; ii<tab->n_a; ii++){

			rownum = ii+1; // number of the row we read in the fits table

			assert(tab->dat[ii]==NULL);
			tab->dat[ii] = load_single_lpDat(fptr,tab->n_h, tab->n_rad, rownum, status);
			if (*status!=EXIT_SUCCESS){
				RELXILL_ERROR("failed to load data from the lp table into memory",status);
				break;
			}
		}

	} while(0);

	if (*status==EXIT_SUCCESS){
		// assigne the value
		(*inp_tab) = tab;
	} else {
		free_lpTable(tab);
	}
	free(fullfilename);

	if (fptr!=NULL) {fits_close_file(fptr,status);}

	return;
}

static void free_relDat(relDat* dat, int nr){
	if (dat!=NULL){
		int ii;
		for (ii=0; ii<nr; ii++){
			if (dat->cosne1 !=NULL) free(dat->cosne1[ii]);
			if (dat->cosne2 !=NULL) free(dat->cosne2[ii]);
			if (dat->trff1 !=NULL) free(dat->trff1[ii]);
			if (dat->trff2 !=NULL) free(dat->trff2[ii]);
		}
		free(dat->cosne1);
		free(dat->cosne2);
		free(dat->trff1);
		free(dat->trff2);

		free(dat->r);
		free(dat->gmin);
		free(dat->gmax);
	}
}

void free_relTable(relTable* tab){
	if(tab!=NULL){
		if (tab->arr!=NULL){
			int ii;
			for (ii=0; ii<tab->n_a; ii++){
				if (tab->arr[ii] !=NULL){
					int jj;
					for (jj=0; jj<tab->n_mu0; jj++){
						free_relDat(tab->arr[ii][jj],tab->n_r);
						free(tab->arr[ii][jj]);
					}
					free(tab->arr[ii]);
				}
			}
			free(tab->arr);
		}
		free(tab->a);
		free(tab->mu0);
		free(tab);
	}
}

lpTable* new_lpTable(int n_a,int n_h, int n_rad, int* status){
	lpTable* tab = (lpTable*) malloc (sizeof(lpTable));
	CHECK_MALLOC_RET_STATUS(tab,status,NULL);

	tab->n_a = n_a;
	tab->n_h = n_h;
	tab->n_rad = n_rad;

	tab->a = NULL;

	tab->dat = (lpDat**) malloc (sizeof(lpDat*)*tab->n_a);
	CHECK_MALLOC_RET_STATUS(tab->dat,status,tab);

	int ii;
	for (ii=0; ii<tab->n_a; ii++){
		tab->dat[ii] = NULL;
	}
	return tab;
}

/* destroy the LP table structure */
void free_lpDat(lpDat* dat,int nh){
	if (dat!=NULL){
		int ii;
		for (ii=0;ii<nh;ii++){
			if (dat->del!=NULL) free(dat->del[ii]);
			if (dat->del_inc!=NULL) free(dat->del_inc[ii]);
			if (dat->intens!=NULL) free(dat->intens[ii]);
		}
		free(dat->del);
		free(dat->del_inc);
		free(dat->intens);

		free(dat->h);
		free(dat->rad);

	}
}

/* destroy the LP table structure */
void free_lpTable(lpTable* tab){
	if (tab!=NULL){
		if (tab->dat!=NULL){
			int ii;
			for (ii=0;ii<tab->n_a;ii++){
				free_lpDat(tab->dat[ii],tab->n_h);
				free(tab->dat[ii]);
			}
			free(tab->dat);
		}
		free(tab->a);
		free(tab);
	}
}

