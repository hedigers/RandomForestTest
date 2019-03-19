
me_test <- function(X,Y,alpha=0.05){
x <<- X
y <<- Y
alphame<<-alpha
py_run_string("
import numpy as np
import freqopttest.util as util
import freqopttest.data as data
import freqopttest.kernel as kernel
import freqopttest.tst as tst
import freqopttest.glo as glo
import sys
#sample source
#n = 3000
#dim = 10
seed = 17
#ss = data.SSGaussMeanDiff(dim, my=1)
#ss = data.SSGaussVarDiff(dim)
#ss = data.SSSameGauss(dim)
#ss = data.SSBlobs()
#dim = ss.dim()
#tst_data = ss.sample(n, seed=seed)
tst_data = data.TSTData(r.x,r.y)
tr, te = tst_data.split_tr_te(tr_proportion=0.5, seed=10)

J = 5
alpha = r.alphame
# seed = 17
op = {'n_test_locs': J, 'seed': seed, 'max_iter': 200, 'locs_step_size': 5.0,
'gwidth_step_size': 0.2, 'tol_fun': 1e-3}
# optimize on the training set
test_locs, gwidth, info = tst.MeanEmbeddingTest.optimize_locs_width(tr, alpha, **op)

met_opt = tst.MeanEmbeddingTest(test_locs, gwidth, alpha)
output = met_opt.perform_test(te)")

decision <- py$output$h0_rejected

return(as.numeric(decision))
}



## Parameters taken from
# https://github.com/wittawatj/interpretable-test/blob/master/freqopttest/ex/ex1_power_vs_n.py
# job_met_opt10




