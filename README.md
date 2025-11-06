# Probabilistic embedding, clustering, and alignment for integrating spatial transcriptomics data with PRECAST

## Simulation 
The simulated dataset examples were in ./Simulation/simulated_datasets folder.

Brief descriptions of simulated scripts (./Simulation/Rcode folder):


**simulated_Scenario1_low(middle,high).R**: Integration data analysis for Scenario 1 with three different batch effects' scales.


**simulated_Scenario2_low(middle,high).R**: Integration data analysis for Scenario 2 with three different batch effects' scales.

**simulated_Scenario3.R**: Integrative data analysis for Scenario 3.

**simulated_Scenario4.R**: Integrative data analysis for Scenario 4.

**simulated_Scenario5.R**: Integrative data analysis for Scenario 5.

## Real data analysis


Brief descriptions of real data analysis scripts (Real_data_analysis folder):

**HCC and mouse liver data**: Four HCC sections are available in the `data` folder. The count matrix is stored in the SeuratObject `HCCx_seu.RDS`, while the metadata and manual annotations (`manual.annotation`) are saved in `meta.data_HCCx.csv`. Additionally, eight mouse liver sections are stored in the SeuratObject `seulist_mouseLiverST8.RDS`, which can be directly loaded into R using the `readRDS()` function.


**dorsolateral_prefrontal_cortex.R**: Integration data analysis for  human dorsolateral prefrontal cortex Visium data

**mouse_liver.R**: Integration data analysis for  mouse liver ST data

**olfactory_bulb.R**:   Integration data analysis for   mouse olfactory bulb Slide-seqV2 data

**hepatocellular_carcinoma.R**: Integration data analysis for  hepatocellular carcinoma Visium data



## Real data results 
The real data results  were visualized with ggplot2 package (Real_data_results folder).

Brief descriptions of real data viualization  scripts:

**dorsolateral_prefrontal_cortex_viualize.R**: Visualization  for  human dorsolateral prefrontal cortex Visium data

**mouse_liver_viualize.R**: Visualization  for  mouse liver ST data

**olfactory_bulb_viualize.R**:  Visualization  for   mouse olfactory bulb Slide-seqV2 data

**hepatocellular_carcinoma_viualize.R**: Visualization  for  hepatocellular carcinoma Visium data
