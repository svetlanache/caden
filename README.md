# caden
 
A package for cross-validated adaptive enrichment signature design (CVASD).<br/>
The package provides three main funcitons:
‘simulate_data’ ‘analyse_simdata’ and ‘analyse_realdata’.
Additional functions are ‘get_stage2’ for simulating a data set for stage 2 based on a data set from stage 2. <br/>
simulate_data:<br/>
     The function simulates covariates data and binary responses for stage 2 to be
     used in the analysis of the CADEN design.<br/>
analyse_simdata:<br/>
     The function computes the power of the design for the simulated
     data according to the CADEN design.<br/>
analyse_realdata:<br/>
     The function analyses the real data according to the CADEN design.<br/>

# istallation
library('devtools')<br/>
install_github('svetlanache/caden')
