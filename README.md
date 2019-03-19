
# Two Sample Testing Based on Random Forest (Replication Code)

We follow the line of using classifiers for two-sample testing and propose several tests based on the Random Forest classifier. The developed tests are easy to use, require no tuning and are applicable for \emph{any} distribution on $\R^p$, even in high-dimensions. We provide a comprehensive treatment for the use of classification for two-sample testing, derive the distribution of our tests under the Null and provide a power analysis, both in theory and with simulations.

Link to the pre-print:
https://arxiv.org/abs/1903.06287

## How to launch the simulation files 

- Install the Python module found on https://github.com/wittawatj/interpretable-test from Wittawat Jitkrittum via the command:
''pip install git+https://github.com/wittawatj/interpretable-test''
Once installed, you should be able to ''import freqopttest''.
- Install our R-package hypoRF found in \hypoRF_Code\hypoRF via R command ''install.packages("../hypoRF/", repos = NULL, type = "source")''.
- Make sure you have installed the R-packages ''ranger'', ''mvtnorm'', ''reticulate'', ''MASS'' and ''parSim'' (only necessary if you are interested in parallel computing).
- Launch the desired simulation R-file locally in the folder \hypoRF_Code\Replicate_Simulations.
If you want to save or plot any result, use the outcommented code at the end of each script.
