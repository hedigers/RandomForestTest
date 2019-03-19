
mmd_test <- function(X,Y, alpha=0.05){
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
#sample source 
#n = 3000
#dim = 10
#seed = 17
alpha=r.alphammd
tst_data = data.TSTData(r.x,r.y)
tr, te = tst_data.split_tr_te(tr_proportion=1, seed=10)
med = util.meddistance(tr.stack_xy(), 1000)
k = kernel.KGauss(med)
mmd_test = tst.QuadMMDTest(k, n_permute=400, alpha=alpha)
output = mmd_test.perform_test(tr)")
                
  
  decision <- py$output$h0_rejected
  
  return(as.numeric(decision))
}





