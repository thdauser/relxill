# -*- mode: Make -*-

refdir_lmod = refdata_localModels
refdir_xilltab = refdata_xillverTable

#### testing ####
test:
	make test-refdata
	make test-functionality

test-refdata:
	make clean-test-refdata
	./test_refdata_relxill.sl 

test-functionality:
	./test_model_functionality.sl

#### clean output data ### 
clean-test-refdata:
	rm -f refdata_*/*/*.pdf
	rm -f refdata_*/*/*.par
	rm -f refdata_*/*/*.tmp

####  Create the reference data #####
refdata:
	make refdata-lmod-create
	make refdata-xilltable-create

refdata-lmod:
	make refdata-lmod-clean
	./create_reference_data.sl

refdata-xilltable:
	make refdata-xilltab-clean
	rm $(refdir_xilltab)/*.fits
	./create_reference_data.sl


refdata-lmod-clean:
	rm -f $(refdir_lmod)/*.fits

refdata-xilltab-clean:
	rm -f $(refdir_xilltab)/*.fits
