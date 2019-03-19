
##############################################################################################################
################################### HOW TO LAUNCH THE REPLICATION R-FILES ####################################
##############################################################################################################

1.) Install the Python module found on https://github.com/wittawatj/interpretable-test from Wittawat Jitkrittum via the command:
	pip install git+https://github.com/wittawatj/interpretable-test
    Once installed, you should be able to import freqopttest.
2.) Install our R-package hypoRF found in \hypoRF_Code\hypoRF via R command install.packages("../hypoRF/", repos = NULL, type = "source").
3.) Make sure you have installed the R-packages ranger, mvtnorm, reticulate, MASS and parSim (only necessary if you are interested in parallel computing).
4.) Launch the desired simulation R-file locally in the folder \hypoRF_Code\Replicate_Simulations.
    If you want to save or plot any result, use the outcommented code at the end of each script.