# Learn2Evaluate
Contains the main function Learn2Evaluate (depends on side functions Learn2Evaluate.R). 
The function Learn2Evaluate allows users to estimate the predictive performance, 
complemented by a confidence bound, by using learning curves. 
For now, this is limited to the AUC as performance measure, and one may only apply ridge regression, lasso regression, or random forest. 

The respository is organized as follows:
1. R scripts: contains the main function Learn2Evaluate, the side functions on which the main function relies, and an application to RNAseq data (run_Learn2Evaluate_real_data.R)
2. Data: contains the RNAseq data, also used in the manuscript (application section)
3. Simulation Results: contains the scripts used to obtain simulation results. These results include the point estimates, confidence bounds, and true AUC values.
   These results are then used to determine the quality of point estimates (quantified by RMSE) and the coverage of the confidence bounds.
