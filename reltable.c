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

#include "reltable.h"

static relDat* new_relDat(int nr, int ng, int* status){
	relDat* dat = (relDat*) malloc (sizeof(relDat));
	CHECK_MALLOC_RET_STATUS(dat,status,dat);

	dat->r = (float*) malloc( sizeof(float) * nr);
	CHECK_MALLOC_RET_STATUS(dat->r,status,dat);
	dat->gmin = (float*) malloc( sizeof(float) * nr);
	CHECK_MALLOC_RET_STATUS(dat->r,status,dat);
	dat->gmax = (float*) malloc( sizeof(float) * nr);
	CHECK_MALLOC_RET_STATUS(dat->r,status,dat);

	dat->cosne1 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->cosne1,status,dat);
	dat->cosne2 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->cosne2,status,dat);
	dat->trff1 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->trff1,status,dat);
	dat->trff2 = (float**) malloc( sizeof(float*) * nr);
	CHECK_MALLOC_RET_STATUS(dat->trff2,status,dat);


	for (int ii=0; ii<nr; ii++){
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

	for (int ii=0; ii<tab->n_a; ii++){

		tab->arr[ii] = (relDat**) malloc (sizeof(relDat*)*tab->n_mu0);
		CHECK_MALLOC_RET_STATUS(tab->arr[ii],status,tab);

		for (int jj=0; jj<tab->n_mu0; jj++){
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

    for (int ii=0; ii<n1;ii++){
        if(fits_read_col(fptr, TFLOAT, colnum, ii+1, 1, nelem ,&nullval,val[ii], &anynul, status)) return;

    }
}

/** load one single data extension from the relline table   */
static relDat* load_single_relDat(fitsfile* fptr, char* extname, int* status){

	int extver = 0;
	fits_movnam_hdu(fptr, BINARY_TBL, extname, extver ,status);
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

		// should be set by previou s routine
		assert(tab!=NULL);
		assert(tab->arr!=NULL);

		// get the full filename
		if (asprintf(&fullfilename, "%s/%s", RELXILL_TABLE_PATH,filename) == -1){
			RELXILL_ERROR("failed to construct full path the rel table",status);
			break;
		}

		// open the file
		if (fits_open_table(&fptr, fullfilename, READONLY, status)) {
			CHECK_RELXILL_ERROR("opening of the rel table failed",status);
			printf("    full path given: %s \n",fullfilename);
			break;
		}

		// first read the axes of the table
		get_reltable_axis(tab->n_a, &(tab->a), "a", "a", fptr, status);
		CHECK_RELXILL_ERROR("reading of spin axis failed",status);

		get_reltable_axis(tab->n_mu0, &(tab->mu0), "mu0", "mu0", fptr, status);
		CHECK_RELXILL_ERROR("reading of mu0 axis failed",status);

		//now load the full table (need to go through all extensions)
		for (int ii=0; ii<tab->n_a; ii++){
			for (int jj=0; jj<tab->n_mu0; jj++){

				if (asprintf(&extname, "%i_%i", ii+1,jj+1) == -1){
					RELXILL_ERROR("failed to construct full path the rel table",status);
					break;
				}

				assert(tab->arr[ii][jj]==NULL);
				tab->arr[ii][jj] = load_single_relDat(fptr, extname, status);
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



static void free_relDat(relDat* dat, int nr){
	if (dat!=NULL){
		for (int ii=0; ii<nr; ii++){
			free(dat->cosne1[ii]);
			free(dat->cosne2[ii]);
			free(dat->trff1[ii]);
			free(dat->trff2[ii]);
		}
		free(dat->cosne1);
		free(dat->cosne2);
		free(dat->trff1);
		free(dat->trff2);

		free(dat->r);
		free(dat->gmin);
		free(dat->gmax);
		free(dat);
	}
}

void free_relTable(relTable* tab){
	if(tab!=NULL){
		if (tab->arr!=NULL){
			for (int ii=0; ii<tab->n_a; ii++){
				if (tab->arr[ii] !=NULL){
					for (int jj=0; jj<tab->n_mu0; jj++){
						free_relDat(tab->arr[ii][jj],tab->n_r);
					}
					free(tab->arr[ii]);
				}
			}
			free(tab->arr);
		}
		free(tab->a);
		free(tab->mu0);
	}
	free(tab);
}
