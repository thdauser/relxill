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
#include "relbase.h"
#include "relutility.h"
#include "reltable.h"
#include "test_relxill.h"
#include "test_xilltab.c"
#include "test_rellp.h"
#include "test_std_functions.h"

void printVersionNumber(int *status) {
  char* buf;
  get_version_number(&buf, status);
  printf("%s",buf);
  free(buf);
}

int main(int argc, char *argv[]){

	int status = EXIT_SUCCESS;



	int do_all = 1;
	int do_relline = 0;
	int do_rellinelp = 0;
	int do_relxill = 0;
	int do_relxilllp = 0;
	int do_relxilllpion = 0;
	int do_relxilldens = 0;
	int do_relxillns = 0;
    int do_relxillco = 0;
	int do_relxilllpdens = 0;
	int do_relxillnthcomp= 0;
	int do_relxilllpnthcomp = 0;
	int do_relxilllpionnthcomp = 0;
	int do_relconv = 0;
	int do_xillver = 0;
	int do_emisTest = 0;

	if (argc>=2){
		if (strcmp(argv[1],"version")==0){
          printVersionNumber(&status);
          return status;
		}

        if (strcmp(argv[1], "relline") == 0) {
            do_relline = 1;
            do_all = 0;
		} else if (strcmp(argv[1],"relconv")==0){
				do_relconv=1;
				do_all=0;
		} else if (strcmp(argv[1],"rellinelp")==0){
				do_rellinelp=1;
				do_all=0;
        } else if (strcmp(argv[1], "relxill") == 0) {
            do_relxill = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilllp") == 0) {
            do_relxilllp = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilldens") == 0) {
            do_relxilldens = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxillNS") == 0) {
            do_relxillns = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxillCO") == 0) {
            do_relxillco = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilllpdens") == 0) {
            do_relxilllpdens = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxillCp") == 0) {
            do_relxillnthcomp = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilllpCp") == 0) {
            do_relxilllpnthcomp = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilllpion") == 0) {
            do_relxilllpion = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "relxilllpionCp") == 0) {
            do_relxilllpionnthcomp = 1;
            do_all = 0;
        } else if (strcmp(argv[1], "xillver") == 0) {
          do_xillver = 1;
          do_all = 0;
        } else if (strcmp(argv[1], "emisTest") == 0) {
          do_emisTest = 1;
          do_all = 0;
        }


    }

	int n = 1;
	if (argc==3){
        n = (int) strtod(argv[2], NULL);
	}

	do{
      print_version_number(&status);

      if (do_all) {
        testStdFunctions(&status);
      }

      if (do_all || do_emisTest){
          test_rellp(&status);
        CHECK_STATUS_BREAK(status);
      }

      if (do_all){


        test_xilltables();
      }


      if (do_all || do_relline){
        status=EXIT_SUCCESS;
        std_eval_relline(&status,n);
        if (status==EXIT_SUCCESS) {
          printf("     ---> successful \n");
        }
      }

      if (do_all || do_rellinelp){
        status=EXIT_SUCCESS;
        std_eval_relline_lp(&status,1);
        if (status==EXIT_SUCCESS) {
          printf("     ---> successful \n");
        }
      }
      
        if (do_all || do_relxill){
            status=EXIT_SUCCESS;
            std_eval_relxill(&status,n);
            CHECK_STATUS_BREAK(status)
            if (status==EXIT_SUCCESS) {
                printf("     ---> successful \n");
            }
            bugtest_eval_relxill(&status);
            CHECK_STATUS_BREAK(status)
            if (status==EXIT_SUCCESS) {
                printf("     ---> successful \n");
            }
        }

        if (do_relconv){
			std_eval_relconv(&status,1);
            CHECK_STATUS_BREAK(status)
			std_eval_relconvlp(&status,1);
            CHECK_STATUS_BREAK(status)
			printf("     ---> successful \n");

		}

		if (do_all | do_xillver ) {
          test_xilltables();
          CHECK_STATUS_BREAK(status)
          std_eval_xillver(&status,1);
          CHECK_STATUS_BREAK(status)
          std_eval_xillver_nthcomp(&status,1);
          CHECK_STATUS_BREAK(status)
          std_eval_xillver_dens_nthcomp(&status,1);
          CHECK_STATUS_BREAK(status)
          printf("     ---> successful \n");

		}


		if (do_all || do_relxillns){
			std_eval_relxill_ns(&status,n);
			CHECK_STATUS_BREAK(status);
			printf("     ---> successful \n");
		}

        if (do_all || do_relxillco) {
            std_eval_relxill_co(&status, n);
            CHECK_STATUS_BREAK(status);
            printf("     ---> successful \n");
        }

        if (do_all || do_relxilllp) {
			std_eval_relxilllp(&status,n);
            CHECK_STATUS_BREAK(status)
			printf("     ---> successful \n");
		}

      if (do_all || do_relxillnthcomp){
          std_eval_relxill_nthcomp(&status,n);
          CHECK_STATUS_BREAK(status)
          std_eval_relxilldens_nthcomp(&status,n);
          CHECK_STATUS_BREAK(status)
          printf("     ---> successful \n");
		}

		if (do_all || do_relxilllpnthcomp){
			std_eval_relxilllp_nthcomp(&status,n);
            CHECK_STATUS_BREAK(status)
          std_eval_relxilllpdens_nthcomp(&status,n);
          CHECK_STATUS_BREAK(status)
			printf("     ---> successful \n");
		}

		if (do_all || do_relxilldens){
			std_eval_relxilldens(&status,n);
            CHECK_STATUS_BREAK(status)
			printf("     ---> successful \n");
		}

		if (do_all || do_relxilllpdens){
			std_eval_relxilllpdens(&status,n);
            CHECK_STATUS_BREAK(status)
            printf("     ---> successful \n");
        }

        if (do_all || do_relxilllpion) {
            std_eval_relxilllpion(&status, n);
            CHECK_STATUS_BREAK(status)
            printf("     ---> successful \n");

        }

        if (do_all || do_relxilllpionnthcomp) {
            std_eval_relxilllpion_nthcomp(&status, n);
            CHECK_STATUS_BREAK(status)
            printf("     ---> successful \n");

        }


		free_cached_tables();
        free_cache();

	} while(0);

	if(status!=EXIT_SUCCESS){
		printf("\n### TESTING NOT SUCCESSFUL \n");
		// free tables
		free_cached_tables();
		free_cache();
	} else {
      printf("\n### TESTING SUCCESSFUL \n");
	}

  return status;
}
