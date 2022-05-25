! CICERO SCM version 9.0 
! Model run '../masan_mess/output_test_rf/ssp245_conc' at 2022.05.25 11:41:34  
output_prefix '../masan_mess/output_test_rf/ssp245_conc'
gaspam_file '../input_RCP/gases_v1RCMIP.txt'
concentration_file 'input/ssp245_conc_RCMIP.txt'
emission_file 'input/ssp245_em_RCMIP.txt'
scenario_file 'input/ssp245_em_RCMIP.txt'
model_end 2100
emission_start 1850
scenario_start 2015
scenario_end 2100
lambda  0.54
bmb_forc -0.01
tropo3_forc  0.40
dirso2_forc -0.46
indso2_forc -0.51
bc_forc  0.20
oc_forc -0.10
ch4_lifetime_mode 'TAR'
model_end 2100
post_scen_mode 'ANNUAL'
post_scen ALL  0.00
