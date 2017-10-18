# Insight_Microbiome_Simulation
This code is meant to recreate some of the work done in the paper "INFANT WEIGHT GAIN TRAJECTORIES LINKED TO ORAL MICROBIOME COMPOSITION" by Craig et. al.  This work is maintained by Dr. Matthew Reimherr at Penn State, all questions about the code or algorithms should be sent to him (mreimherr@psu.edu).  

### Files
The main codes are "ex1_combine_abundance_CS.R" and "ex_flame_CS_and_reg.R".  The first goes through our algorithm for merging bacterial abundances in the gut microbiome of our subjects (in this case the children being monitored).  The second carries out the variable selection using FLAME (a LASSO style procedure that allows the outcomes to be functions) as well as a functional regression example which creates an artificial signal for the curves and then estimates it.  The necessary R scripts which are used by the primary code are in the "RFunctions" folder.  The data for the functiona regressions is stored in "flame_CS_sim_data.RData", while the data and output for the abundances is in the "data" folder 