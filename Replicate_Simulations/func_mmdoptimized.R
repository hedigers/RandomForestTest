
mmdoptimized_test <- function(X,Y, alpha=0.05){
x <<- X
y <<- Y
alphammd <<- alpha
py_run_string("
import numpy as np
import matplotlib.pyplot as plt
import freqopttest.util as util
import freqopttest.data as data
import freqopttest.kernel as kernel
import freqopttest.tst as tst
import freqopttest.glo as glo
import scipy.stats as stats
import sys
alpha=r.alphammd
tst_data = data.TSTData(r.x,r.y)
tr, te = tst_data.split_tr_te(tr_proportion=0.5, seed=10)
med = util.meddistance(tr.stack_xy(), 1000)
list_gwidth = np.hstack( ( (med**2) *(2.0**np.linspace(-4, 4, 30) ) ) )
list_gwidth.sort()
list_kernels = [kernel.KGauss(gw2) for gw2 in list_gwidth]
besti, powers = tst.QuadMMDTest.grid_search_kernel(tr, list_kernels, alpha)
best_ker = list_kernels[besti]
mmd_test = tst.QuadMMDTest(best_ker, n_permute=400, alpha=alpha)
output = mmd_test.perform_test(te)")
  
  
  decision <- py$output$h0_rejected
  
  return(as.numeric(decision))
}

## Parameters taken from
# https://github.com/wittawatj/interpretable-test/blob/master/freqopttest/ex/ex1_power_vs_n.py
# job_quad_mmd






